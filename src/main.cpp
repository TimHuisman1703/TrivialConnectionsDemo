// Disable compiler warnings in third-party code (which we cannot change).
#include <framework/disable_all_warnings.h>
#include <framework/opengl_includes.h>
DISABLE_WARNINGS_PUSH()
#include <fmt/format.h>
// Include glad before glfw3
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/rotate_vector.hpp>
#include <glm/mat4x4.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
// #define STB_IMAGE_IMPLEMENTATION
#include <stb/stb_image.h>
DISABLE_WARNINGS_POP()
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib> // EXIT_FAILURE
#include <framework/mesh.h>
#include <framework/shader.h>
#include <framework/trackball.h>
#include <framework/window.h>
#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw.h>
#include <imgui/imgui_impl_opengl2.h>
#include <iostream>
#include <numeric>
#include <optional>
#include <span>
#include <stack>
#include <vector>

#include "buffer.h"
#include "mesh_ex.h"
#include "solve.h"
#include "tree_cotree.h"

// Configuration
const int WIDTH = 1200;
const int HEIGHT = 800;

enum Stage {
	CHOOSE_MESH,
	VERTEX_K,
	TREE_COTREE,
	INSPECT_RESULT,

	STAGE_NR_ITEMS
};

// Copied (and adapted) from 3DCGA Exercise Set 03
static std::optional<glm::vec3> getWorldPositionOfPixel(const Trackball& trackball, const glm::vec2& pixel)
{
	float depth;
	glReadPixels(static_cast<int>(pixel.x), static_cast<int>(pixel.y), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &depth);

	if (depth == 1.0f) {
		// This is a work around for a bug in GCC:
		// https://gcc.gnu.org/bugzilla/show_bug.cgi?id=80635
		//
		// This bug will emit a warning about a maybe uninitialized value when writing:
		// return {};
		constexpr std::optional<glm::vec3> tmp;
		return tmp;
	}

	// Coordinates convert from pixel space to OpenGL screen space (range from -1 to +1)
	const glm::vec3 win{ pixel, depth };

	// View matrix
	const glm::mat4 view = trackball.viewMatrix();
	const glm::mat4 projection = trackball.projectionMatrix();

	const glm::vec4 viewport{ 0, 0, WIDTH, HEIGHT };
	return glm::unProject(win, view, projection, viewport);
}
static size_t getClosestVertexIndex(const MeshEx& mesh_ex, const glm::vec3& pos)
{
	const auto iter = std::min_element(
		std::begin(mesh_ex.vertices), std::end(mesh_ex.vertices),
		[&](const VertexEx& lhs, const VertexEx& rhs) {
			return glm::length(lhs.position - pos) < glm::length(rhs.position - pos);
		});
	return (size_t)std::distance(std::begin(mesh_ex.vertices), iter);
}
static size_t getClosestEdgeIndex(const MeshEx& mesh_ex, const glm::vec3& pos)
{
	const auto iter = std::min_element(
		std::begin(mesh_ex.edges), std::end(mesh_ex.edges),
		[&](const EdgeEx& lhs, const EdgeEx& rhs) {
			glm::vec3 lhs_v_a_pos = mesh_ex.vertices[lhs.vertices[0]].position;
			glm::vec3 lhs_v_b_pos = mesh_ex.vertices[lhs.vertices[1]].position;
			glm::vec3 rhs_v_a_pos = mesh_ex.vertices[rhs.vertices[0]].position;
			glm::vec3 rhs_v_b_pos = mesh_ex.vertices[rhs.vertices[1]].position;
			return glm::length((lhs_v_a_pos + lhs_v_b_pos) * 0.5f - pos) < glm::length((rhs_v_a_pos + rhs_v_b_pos) * 0.5f - pos);
		});
	return (size_t)std::distance(std::begin(mesh_ex.edges), iter);
}

// Program entry point. Everything starts here.
int main(int argc, char** argv)
{
	// Window settings
	Window window{ "Trivial Connections - Demonstration", glm::ivec2(WIDTH, HEIGHT), OpenGLVersion::GL45 };

	// Camera settings
	glm::vec3 look_at{ 0.0f };
	glm::vec3 rotations{ 0.6f, 0.3f, 0.0f };
	float fov_y = 30.0f;
	float dist = 4.5f;
	Trackball trackball{ &window, glm::radians(fov_y) };
	trackball.setCamera(look_at, rotations, dist);

	// Create buffers
	MeshBuffer mesh_buffer_dual_edges = MeshBuffer::empty();
	MeshBuffer mesh_buffer_tree = MeshBuffer::empty();
	MeshBuffer mesh_buffer_cotree = MeshBuffer::empty();
	std::vector<MeshBuffer> mesh_buffers_noncon;
	MeshBuffer mesh_buffer_field = MeshBuffer::empty();
	MeshBuffer mesh_buffer_travel_vector_original = MeshBuffer::empty();
	MeshBuffer mesh_buffer_travel_field_original = MeshBuffer::empty();
	MeshBuffer mesh_buffer_travel_vector_adjusted = MeshBuffer::empty();
	MeshBuffer mesh_buffer_travel_field_adjusted = MeshBuffer::empty();

	// Build shaders
	const Shader normal_shader = ShaderBuilder()
		.addStage(GL_VERTEX_SHADER, RESOURCE_ROOT "shaders/normal.vert")
		.addStage(GL_FRAGMENT_SHADER, RESOURCE_ROOT "shaders/normal.frag")
		.build();
	const Shader wireframe_shader = ShaderBuilder()
		.addStage(GL_VERTEX_SHADER, RESOURCE_ROOT "shaders/wireframe.vert")
		.addStage(GL_FRAGMENT_SHADER, RESOURCE_ROOT "shaders/wireframe.frag")
		.build();
	const Shader singularity_shader = ShaderBuilder()
		.addStage(GL_VERTEX_SHADER, RESOURCE_ROOT "shaders/singularity.vert")
		.addStage(GL_FRAGMENT_SHADER, RESOURCE_ROOT "shaders/singularity.frag")
		.build();

	// Mesh data
	MeshBuffer mesh_buffer_main;
	MeshBuffer mesh_buffer_wireframe;
	MeshEx mesh_ex;

	// Load singularity mesh
	const Mesh mesh_singularity = loadMesh(std::string(RESOURCE_ROOT) + "resources/singularity.obj")[0];
	const MeshBuffer mesh_buffer_singularity = MeshBuffer::fromTriangles(mesh_singularity.vertices, mesh_singularity.triangles);

	// Parameters
	bool show_imgui = true;
	Stage stage = CHOOSE_MESH;
	int selected_vertex_idx = 0;
	int clicking_edge = 0;
	bool show_tree = true;
	bool show_cotree = true;
	bool show_noncons = true;
	int selected_noncon_item_idx = 0;
	std::vector<std::string> noncon_items;
	float travel_speed = 0.01f;
	float start_angle = 0.0f;
	bool show_travel_original = true;
	bool show_travel_adjusted = true;
	bool show_singularity_data = true;

	// Algorithm data
	std::vector<int> tree_assignment;
	std::vector<std::vector<int>> noncon_cycles;
	std::vector<int> noncon_ks;
	std::vector<double> adjustment_angles;

	// Cycle data
	std::vector<bool> dual_edges_selected;
	bool dual_edges_reload = true;
	std::vector<int> dual_edges_path;
	bool is_noncon;
	double curvature_partition;
	double curvature_partition_adjusted;

	// Field
	std::vector<double> field_angles;

	// Travel
	std::vector<glm::vec3> travel_path;
	std::vector<glm::vec3> travel_normals;
	std::vector<double> travel_angles_original;
	std::vector<double> travel_angles_adjusted;
	std::vector<double> travel_angle_adjustments;
	int travel_index;
	float travel_frac;

	// Handle key press
	window.registerKeyCallback([&](int key, int scancode, int action, int mods) {
		if (key == '\\' && action == GLFW_PRESS) {
			show_imgui = !show_imgui;
		}

		if (key == GLFW_KEY_LEFT_SHIFT && action == GLFW_RELEASE)
			clicking_edge = 0;

		if (action != GLFW_RELEASE)
			return;
		});

	// Handle mouse 
	window.registerMouseButtonCallback([&](int button, int action, int mods) {
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE)
			clicking_edge = 0;

		if (!window.isKeyPressed(GLFW_KEY_LEFT_SHIFT))
			return;

		if (stage == VERTEX_K) {
			if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
				const auto opt_world_point = getWorldPositionOfPixel(trackball, window.getCursorPixel());
				if (opt_world_point)
					selected_vertex_idx = getClosestVertexIndex(mesh_ex, *opt_world_point);
				return;
			}
		}
		else if (stage == INSPECT_RESULT) {
			if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
				clicking_edge = 2;
				return;
			}
		}
		});

	glEnable(GL_DEPTH_TEST);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_BLEND);

	const std::array<std::string, 4> stage_names{
		"Choosing Mesh",
		"Setting Singularities",
		"Tree-Cotree Decomposition",
		"Inspect Results",
	};
	const std::array<std::string, 4> stage_descriptions{
		"Select a mesh to construct a\n  vector field on.",
		"Set the singularities of each\n  vertex.\nSelect a vertex by holding\n  SHIFT and clicking on the\n  mesh.\nIncrease or decrease its index\n  in the GUI.",
		"View the tree and cotree\n  generated for this mesh.\nHighlight a noncontractible\n  cycle in the GUI.\nIncrease or decrease its index\n  in the GUI.",
		"Inspect the vector field.\nDraw a cycle by selecting a\n  set of edges.\nSelect an edge by holding SHIFT\n  and clicking on the mesh.\nInspect data about the cycle\n  in the GUI.",
	};

	// Method for loading new mesh
	auto loadMainMesh = [&](std::string name) {
		// Load main mesh
		const Mesh mesh = loadMesh(std::string(RESOURCE_ROOT) + fmt::format("resources/{}.obj", name.c_str()))[0];
		mesh_buffer_main = MeshBuffer::fromTriangles(mesh.vertices, mesh.triangles);

		// Create extended mesh representation
		mesh_ex = MeshEx::fromMesh(mesh);
		std::vector<glm::uvec2> mesh_edges = {};
		for (const EdgeEx& e : mesh_ex.edges) {
			mesh_edges.push_back(e.vertices);
		}

		// Create wireframe buffer
		mesh_buffer_wireframe = MeshBuffer::fromEdges(mesh.vertices, mesh_edges);

		tree_assignment = treeCotreeDecompose(mesh_ex);
		noncon_cycles = findNoncontractibleCycles(mesh_ex, tree_assignment);

		// Create tree data for rendering
		std::vector<glm::vec3> tree_positions;
		std::vector<int> tree_indices;
		for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++) {
			if (tree_assignment[e_idx] == 1) {
				const EdgeEx& e = mesh_ex.edges[e_idx];

				glm::vec3 v_a_pos = mesh_ex.vertices[e.vertices[0]].position;
				glm::vec3 v_b_pos = mesh_ex.vertices[e.vertices[1]].position;

				tree_positions.push_back(v_a_pos);
				tree_positions.push_back(v_b_pos);

				tree_indices.push_back(tree_positions.size() - 2);
				tree_indices.push_back(tree_positions.size() - 1);
			}
		}
		mesh_buffer_tree.load(tree_positions, tree_indices);

		// Create cotree data for rendering
		std::vector<glm::vec3> cotree_positions;
		std::vector<int> cotree_indices;
		for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++) {
			if (tree_assignment[e_idx] == -1) {
				const EdgeEx& e = mesh_ex.edges[e_idx];

				glm::vec3 e_center = 0.5f * (mesh_ex.vertices[e.vertices[0]].position + mesh_ex.vertices[e.vertices[1]].position);
				glm::vec3 f_a_center = mesh_ex.centerOfMass(e.faces[0]);
				glm::vec3 f_b_center = mesh_ex.centerOfMass(e.faces[1]);

				cotree_positions.push_back(f_a_center);
				cotree_positions.push_back(e_center);
				cotree_positions.push_back(f_b_center);

				cotree_indices.push_back(cotree_positions.size() - 3);
				cotree_indices.push_back(cotree_positions.size() - 2);
				cotree_indices.push_back(cotree_positions.size() - 2);
				cotree_indices.push_back(cotree_positions.size() - 1);
			}
		}
		mesh_buffer_cotree.load(cotree_positions, cotree_indices);

		// Create noncon-cycle data for rendering
		for (const MeshBuffer& mb : mesh_buffers_noncon)
			mb.cleanUp();
		mesh_buffers_noncon = {};

		for (std::vector<int>& noncon_cycle : noncon_cycles) {
			MeshBuffer mesh_buffer_noncon = MeshBuffer::empty();
			std::vector<glm::vec3> noncon_positions;
			std::vector<int> noncon_indices;

			for (int e_idx : noncon_cycle) {
				const EdgeEx& e = mesh_ex.edges[e_idx];

				glm::vec3 e_center = 0.5f * (mesh_ex.vertices[e.vertices[0]].position + mesh_ex.vertices[e.vertices[1]].position);
				glm::vec3 f_a_center = mesh_ex.centerOfMass(e.faces[0]);
				glm::vec3 f_b_center = mesh_ex.centerOfMass(e.faces[1]);

				noncon_positions.push_back(f_a_center);
				noncon_positions.push_back(e_center);
				noncon_positions.push_back(f_b_center);

				noncon_indices.push_back(noncon_positions.size() - 3);
				noncon_indices.push_back(noncon_positions.size() - 2);
				noncon_indices.push_back(noncon_positions.size() - 2);
				noncon_indices.push_back(noncon_positions.size() - 1);
			}

			mesh_buffer_noncon.load(noncon_positions, noncon_indices);
			mesh_buffers_noncon.push_back(mesh_buffer_noncon);
		}

		// Initialize noncon-cycle k's
		noncon_ks = {};
		for (int i = 0; i < noncon_cycles.size(); i++)
			noncon_ks.push_back(0);

		// Clear selected edges
		dual_edges_selected = {};
		for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++)
			dual_edges_selected.push_back(false);

		// Create GUI selection items
		noncon_items = { "None" };
		for (int i = 0; i < noncon_cycles.size(); i++)
			noncon_items.push_back(std::format("Cycle #{}", i + 1));
	};

	loadMainMesh("icosphere");

	while (!window.shouldClose()) {
		// Update input
		window.updateInput();

		// Poincare-Hopf calculations
		int total_k = 0;
		for (const VertexEx& v : mesh_ex.vertices)
			total_k += v.k;
		for (const int k : noncon_ks)
			total_k += k;
		int euler_mesh_characteristic = mesh_ex.vertices.size() - mesh_ex.edges.size() + mesh_ex.faces.size();

		// Stage checks
		bool prev_enabled = stage > 0;
		bool next_enabled = stage < STAGE_NR_ITEMS - 1;
		if (stage == TREE_COTREE) {
			if (total_k != euler_mesh_characteristic)
				next_enabled = false;
		}

		// Draw GUI
		if (show_imgui) {
			ImGui::Begin("Control");

			float width = ImGui::GetContentRegionAvail().x;
			float halfWidth = (width - ImGui::GetStyle().ItemSpacing.x) / 2;

			ImGui::Text("Stage %i: %s", stage + 1, stage_names[stage].c_str());

			if (ImGui::Button(prev_enabled ? "Previous" : "##", ImVec2(halfWidth, 19))) {
				if (prev_enabled) {
					stage = static_cast<Stage>(stage - 1);

					selected_vertex_idx = -1;
					dual_edges_reload = true;
				}
			}
			ImGui::SameLine();
			if (ImGui::Button(next_enabled ? "Next" : "###", ImVec2(halfWidth, 19))) {
				if (next_enabled) {
					stage = static_cast<Stage>(stage + 1);

					selected_vertex_idx = -1;
					selected_noncon_item_idx = 0;
					dual_edges_reload = true;

					if (stage == INSPECT_RESULT) {
						dual_edges_reload = true;
						std::vector<std::pair<std::vector<int>, int>> cycles = getCycles(mesh_ex, noncon_cycles, noncon_ks);
						adjustment_angles = calculateAdjustmentAngles(mesh_ex, cycles);
					}
				}
			}

			ImGui::Separator();
			ImGui::Text(stage_descriptions[stage].c_str());
			ImGui::Separator();

			if (stage == CHOOSE_MESH) {
				std::string mesh_name;
				if (ImGui::Button("Cube (g = 0)", ImVec2(width, 19)))
					mesh_name = "cube";
				if (ImGui::Button("Icosphere (g = 0)", ImVec2(width, 19)))
					mesh_name = "icosphere";
				if (ImGui::Button("Torus (g = 1)", ImVec2(width, 19)))
					mesh_name = "torus";
				if (ImGui::Button("What? (g = 6)", ImVec2(width, 19)))
					mesh_name = "what";
				if (ImGui::Button("Bunny (g = 0)", ImVec2(width, 19)))
					mesh_name = "bunny";

				if (!mesh_name.empty())
					loadMainMesh(mesh_name);
			}
			else if (stage == VERTEX_K) {
				if (selected_vertex_idx > -1) {
					ImGui::Text("Selected vertex: #%i", selected_vertex_idx);
					ImGui::Text("Singularity index (k):");
					ImGui::Indent();
					ImGui::InputInt("k", &mesh_ex.vertices[selected_vertex_idx].k);
					ImGui::Unindent();
				}
				else {
					ImGui::Text("Selected vertex: None");
				}

				ImGui::Separator();
				if (total_k == euler_mesh_characteristic) {
					ImGui::TextColored(ImVec4(0.2, 1.0, 0.4, 1.0), "Poincare-Hopf satisfied");
					ImGui::Text("  sum(k) = %i == %i = 2-2g", total_k, euler_mesh_characteristic);
				}
				else {
					ImGui::TextColored(ImVec4(1.0, 0.4, 0.2, 1.0), "Poincare-Hopf not satisfied");
					ImGui::Text("  sum(k) = %i != %i = 2-2g", total_k, euler_mesh_characteristic);
					ImGui::TextColored(ImVec4(1.0, 1.0, 0.2, 1.0), "This can still be fixed in the\n  next step");
				}
			}
			else if (stage == TREE_COTREE) {
				ImGui::Checkbox("Show tree", &show_tree);
				ImGui::Checkbox("Show cotree", &show_cotree);
				ImGui::Checkbox("Show noncontractible cycles", &show_noncons);
				if (show_noncons) {
					std::vector<const char*> noncon_items_pointers;
					for (const std::string& item : noncon_items)
						noncon_items_pointers.push_back(item.c_str());
					ImGui::Text("Selected cycle:");
					ImGui::Indent();
					ImGui::Combo("##", &selected_noncon_item_idx, noncon_items_pointers.data(), (int)noncon_items_pointers.size());
					ImGui::Unindent();
					if (selected_noncon_item_idx > 0) {
						ImGui::Text("Singularity index (k):");
						ImGui::Indent();
						ImGui::InputInt("k", &noncon_ks[selected_noncon_item_idx - 1]);
						ImGui::Unindent();
					}
				}
				ImGui::Separator();
				if (total_k == euler_mesh_characteristic) {
					ImGui::TextColored(ImVec4(0.2, 1.0, 0.4, 1.0), "Poincare-Hopf satisfied");
					ImGui::Text("  sum(k) = %i == %i = 2-2g", total_k, euler_mesh_characteristic);
				}
				else {
					ImGui::TextColored(ImVec4(1.0, 0.4, 0.2, 1.0), "Poincare-Hopf not satisfied");
					ImGui::Text("  sum(k) = %i != %i = 2-2g", total_k, euler_mesh_characteristic);
					ImGui::TextColored(ImVec4(1.0, 1.0, 0.2, 1.0), "This must be fixed before\n  proceeding to the next step");
				}
			}
			else if (stage == INSPECT_RESULT) {
				if (ImGui::Button("Clear", ImVec2(width, 19))) {
					dual_edges_selected = {};
					for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++)
						dual_edges_selected.push_back(false);
					dual_edges_reload = true;
				}
				ImGui::SliderFloat("Angle", &start_angle, 0.0f, 2 * glm::pi<float>());
				ImGui::Checkbox("Show singularity data", &show_singularity_data);
				if (!dual_edges_path.empty()) {
					int num_rotations = glm::round(0.5 * curvature_partition_adjusted / glm::pi<double>());

					ImGui::TextColored(ImVec4(0.2, 1.0, 0.4, 1.0), "Cycle complete");
					ImGui::Text(is_noncon ? "  Type       = Noncontractible" : "  Type       = Contractible");
					ImGui::Text("  #Edges     = %i", dual_edges_path.size());
					ImGui::Text("  Rotation");
					ImGui::Text("    Original = %f", curvature_partition);
					ImGui::Text("    Adjusted = %f", curvature_partition_adjusted);
					ImGui::Text("      %i rotation(s)", num_rotations);
					ImGui::Checkbox("Show original vectors", &show_travel_original);
					ImGui::Checkbox("Show adjusted vectors", &show_travel_adjusted);
					ImGui::SliderFloat("Speed", &travel_speed, 0.0f, 0.05f);
				}
				else {
					ImGui::TextColored(ImVec4(1.0, 0.4, 0.2, 1.0), "Cycle incomplete");
				}
			}

			ImGui::End();
			ImGui::Render();
		}

		// Logic
		if (stage == INSPECT_RESULT) {
			if (clicking_edge) {
				const auto opt_world_point = getWorldPositionOfPixel(trackball, window.getCursorPixel());
				if (opt_world_point) {
					int closest_edge_index = getClosestEdgeIndex(mesh_ex, *opt_world_point);

					if (clicking_edge == 2)
						clicking_edge = 1 - 2 * dual_edges_selected[closest_edge_index];

					bool selected = (clicking_edge + 1) / 2;
					if (dual_edges_selected[closest_edge_index] != selected) {
						dual_edges_selected[closest_edge_index] = selected;
						dual_edges_reload = true;
					}
				}
			}

			if (dual_edges_reload) {
				travel_path = {};
				travel_path = {};
				travel_normals = {};
				travel_angles_original = { 0.0f };
				travel_angles_adjusted = { 0.0f };
				travel_angle_adjustments = {};

				// Trace cycle
				dual_edges_path = {};
				std::vector<std::vector<int>> edges_per_face(mesh_ex.faces.size());

				int f_start_idx = -1;
				int e_start_idx = -1;
				int expected_length = 0;
				for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++) {
					if (dual_edges_selected[e_idx]) {
						const EdgeEx& e = mesh_ex.edges[e_idx];
						edges_per_face[e.faces[0]].push_back(e_idx);
						edges_per_face[e.faces[1]].push_back(e_idx);

						if (f_start_idx == -1) {
							f_start_idx = e.faces[0];
							e_start_idx = e_idx;
						}
						expected_length++;
					}
				}

				if (f_start_idx != -1) {
					dual_edges_path = { e_start_idx };

					int f_curr_idx = mesh_ex.otherFace(e_start_idx, f_start_idx);
					int e_curr_idx = e_start_idx;
					while (f_curr_idx != f_start_idx) {
						std::vector<int> edges = edges_per_face[f_curr_idx];

						if (edges.size() != 2) {
							dual_edges_path = {};
							break;
						}

						if (edges[0] != e_curr_idx)
							e_curr_idx = edges[0];
						else
							e_curr_idx = edges[1];

						f_curr_idx = mesh_ex.otherFace(e_curr_idx, f_curr_idx);
						dual_edges_path.push_back(e_curr_idx);
					}

					if (dual_edges_path.size() != expected_length)
						dual_edges_path = {};

					std::vector<bool> dual_edges_vertex_partition = {};
					for (int v_idx = 0; v_idx < mesh_ex.vertices.size(); v_idx++)
						dual_edges_vertex_partition.push_back(false);
					if (!dual_edges_path.empty()) {
						int v_start_idx = mesh_ex.edges[e_curr_idx].vertices[0];

						std::stack<int> stack;
						stack.push(v_start_idx);
						int partition_size = 0;
						while (!stack.empty()) {
							int v_curr_idx = stack.top();
							stack.pop();

							if (dual_edges_vertex_partition[v_curr_idx])
								continue;
							dual_edges_vertex_partition[v_curr_idx] = true;
							partition_size++;

							for (int e_idx : mesh_ex.vertices[v_curr_idx].edges) {
								if (!dual_edges_selected[e_idx])
									stack.push(mesh_ex.otherVertex(e_idx, v_curr_idx));
							}
						}

						is_noncon = partition_size == dual_edges_vertex_partition.size();

						if (2 * partition_size > dual_edges_vertex_partition.size()) {
							for (int v_idx = 0; v_idx < dual_edges_vertex_partition.size(); v_idx++)
								dual_edges_vertex_partition[v_idx] = !dual_edges_vertex_partition[v_idx];
						}

						curvature_partition = mesh_ex.angleOnPath(dual_edges_path);
						curvature_partition_adjusted = mesh_ex.angleOnPathAdjusted(dual_edges_path, adjustment_angles);
					}
				}

				// Create cycle data for rendering
				std::vector<glm::vec3> dual_edges_positions;
				std::vector<int> dual_edges_indices;
				if (dual_edges_path.empty()) {
					// Incomplete path
					for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++) {
						if (dual_edges_selected[e_idx]) {
							const EdgeEx& e = mesh_ex.edges[e_idx];

							glm::vec3 e_center = 0.5f * (mesh_ex.vertices[e.vertices[0]].position + mesh_ex.vertices[e.vertices[1]].position);
							glm::vec3 f_a_center = mesh_ex.centerOfMass(e.faces[0]);
							glm::vec3 f_b_center = mesh_ex.centerOfMass(e.faces[1]);

							const FaceEx& f_a = mesh_ex.faces[e.faces[0]];
							const FaceEx& f_b = mesh_ex.faces[e.faces[1]];

							dual_edges_positions.push_back(f_a_center + glm::normalize(f_a.normal) * 0.005f);
							dual_edges_positions.push_back(e_center + glm::normalize(f_a.normal + f_b.normal) * 0.005f);
							dual_edges_positions.push_back(f_b_center + glm::normalize(f_b.normal) * 0.005f);
							dual_edges_indices.push_back(dual_edges_positions.size() - 3);
							dual_edges_indices.push_back(dual_edges_positions.size() - 2);
							dual_edges_indices.push_back(dual_edges_positions.size() - 2);
							dual_edges_indices.push_back(dual_edges_positions.size() - 1);
						}
					}

					// Field
					std::vector<bool> seen;
					field_angles = {};
					for (int f_idx = 0; f_idx < mesh_ex.faces.size(); f_idx++) {
						seen.push_back(f_idx == 0);
						field_angles.push_back(0.0);
					}
					std::queue<int> queue = {};
					queue.push(0);
					while (!queue.empty()) {
						int f_curr_idx = queue.front();
						queue.pop();
						const FaceEx& f_curr = mesh_ex.faces[f_curr_idx];

						for (int e_b_idx : f_curr.edges) {
							int f_next_idx = mesh_ex.otherFace(e_b_idx, f_curr_idx);
							if (seen[f_next_idx])
								continue;
							seen[f_next_idx] = true;

							const FaceEx& f_next = mesh_ex.faces[f_next_idx];
							int e_a_idx = f_curr.edges[0];
							int e_c_idx = f_next.edges[0];

							const EdgeEx& e_b = mesh_ex.edges[e_b_idx];
							double angle_adjustment = adjustment_angles[e_b_idx];
							if (e_b.faces[0] == f_next_idx)
								angle_adjustment = -angle_adjustment;

							field_angles[f_next_idx] = field_angles[f_curr_idx]
								+ mesh_ex.turnAngleBetweenEdges(e_a_idx, e_b_idx)
								+ glm::pi<float>() - angle_adjustment
								+ mesh_ex.turnAngleBetweenEdges(e_b_idx, e_c_idx);
							queue.push(f_next_idx);
						}
					}
				}
				else {
					// Complete path (tighter)
					for (int i = 0; i < dual_edges_path.size(); i++) {
						int e_a_idx = dual_edges_path[i];
						int e_b_idx = dual_edges_path[(i + 1) % dual_edges_path.size()];
						int e_c_idx = dual_edges_path[(i + 2) % dual_edges_path.size()];
						int f_to_idx = mesh_ex.commonFaceOfEdges(e_a_idx, e_b_idx);
						int f_from_idx = mesh_ex.otherFace(e_a_idx, f_to_idx);
						int f_next_idx = mesh_ex.otherFace(e_b_idx, f_to_idx);

						const EdgeEx& e_a = mesh_ex.edges[e_a_idx];
						const EdgeEx& e_b = mesh_ex.edges[e_b_idx];
						const EdgeEx& e_c = mesh_ex.edges[e_c_idx];
						const FaceEx& f_from = mesh_ex.faces[f_from_idx];
						const FaceEx& f_to = mesh_ex.faces[f_to_idx];
						const FaceEx& f_next = mesh_ex.faces[f_next_idx];

						glm::vec3 e_a_center = 0.5f * (mesh_ex.vertices[e_a.vertices[0]].position + mesh_ex.vertices[e_a.vertices[1]].position);
						glm::vec3 e_b_center = 0.5f * (mesh_ex.vertices[e_b.vertices[0]].position + mesh_ex.vertices[e_b.vertices[1]].position);
						glm::vec3 e_c_center = 0.5f * (mesh_ex.vertices[e_c.vertices[0]].position + mesh_ex.vertices[e_c.vertices[1]].position);

						glm::vec3 dir_ab = glm::normalize(e_b_center - e_a_center);
						glm::vec3 dir_b = glm::normalize(mesh_ex.vertices[e_b.vertices[1]].position - mesh_ex.vertices[e_b.vertices[0]].position);
						if (glm::dot(glm::cross(dir_ab, dir_b), f_to.normal) < 0)
							dir_b = -dir_b;
						glm::vec3 dir_bc = glm::normalize(e_c_center - e_b_center);

						double angle_curve = glm::acos(glm::dot(dir_bc, dir_b)) - glm::acos(glm::dot(dir_ab, dir_b));
						double angle_adjustment = adjustment_angles[e_b_idx];
						if (e_b.faces[0] == f_to_idx)
							angle_adjustment = -angle_adjustment;

						travel_path.push_back(e_a_center);
						travel_normals.push_back(f_to.normal);
						travel_angles_original.push_back(travel_angles_original[travel_angles_original.size() - 1] + angle_curve);
						travel_angles_adjusted.push_back(travel_angles_adjusted[travel_angles_adjusted.size() - 1] + angle_curve + angle_adjustment);
						travel_angle_adjustments.push_back(angle_adjustment);

						dual_edges_positions.push_back(e_a_center + glm::normalize(f_from.normal + f_to.normal) * 0.005f);
						dual_edges_positions.push_back(e_b_center + glm::normalize(f_to.normal + f_next.normal) * 0.005f);
						dual_edges_indices.push_back(dual_edges_positions.size() - 2);
						dual_edges_indices.push_back(dual_edges_positions.size() - 1);
					}
				}
				mesh_buffer_dual_edges.load(dual_edges_positions, dual_edges_indices);

				travel_index = 0;
				travel_frac = 0.5f;

				dual_edges_reload = false;
			}

			if (!travel_path.empty()) {
				// Update travel vector data
				{
					glm::vec3 pos_a = travel_path[travel_index];
					glm::vec3 pos_b = travel_path[(travel_index + 1) % travel_path.size()];
					glm::vec3 normal = travel_normals[travel_index];
					float angle_original = travel_angles_original[travel_index];
					float angle_adjusted = travel_angles_adjusted[travel_index];
					if (travel_index == 0 && travel_frac < 0.5f) {
						angle_original = travel_angles_original[travel_angles_original.size() - 1];
						angle_adjusted = travel_angles_adjusted[travel_angles_adjusted.size() - 1];
					}
					angle_adjusted -= travel_angle_adjustments[(travel_index + travel_angle_adjustments.size() - 1) % travel_angle_adjustments.size()] * glm::max(0.0f, 0.5f - travel_frac);
					angle_adjusted += travel_angle_adjustments[travel_index] * glm::max(0.0f, travel_frac - 0.5f);

					glm::vec3 origin = pos_a * (1.0f - travel_frac) + pos_b * travel_frac + normal * 0.01f;;
					glm::vec3 direction_original = glm::rotate(glm::normalize(pos_b - pos_a), start_angle + angle_original, normal);
					glm::vec3 direction_adjusted = glm::rotate(glm::normalize(pos_b - pos_a), start_angle + angle_adjusted, normal);

					double edge_length = glm::length(pos_b - pos_a);
					travel_frac += travel_speed / edge_length;
					if (travel_frac >= 1.0f) {
						travel_index = (travel_index + 1) % travel_path.size();
						travel_frac = 0.0f;
					}

					// Create travel vector data for rendering
					mesh_buffer_travel_vector_original.load({ origin, origin + direction_original * 0.2f }, { 0, 1 });
					mesh_buffer_travel_vector_adjusted.load({ origin, origin + direction_adjusted * 0.2f }, { 0, 1 });
				}

				// Create travel field data for rendering
				{
					std::vector<glm::vec3> travel_field_original_positions;
					std::vector<int> travel_field_original_indices;
					std::vector<glm::vec3> travel_field_adjusted_positions;
					std::vector<int> travel_field_adjusted_indices;

					for (int i = 0; i <= travel_path.size(); i++) {
						glm::vec3 pos_a = travel_path[i % travel_path.size()];
						glm::vec3 pos_b = travel_path[(i + 1) % travel_path.size()];
						glm::vec3 normal = travel_normals[i % travel_path.size()];
						float angle_original = travel_angles_original[i];
						float angle_adjusted = travel_angles_adjusted[i];

						glm::vec3 origin = 0.5f * (pos_a + pos_b) + normal * 0.01f;
						glm::vec3 direction_original = glm::rotate(glm::normalize(pos_b - pos_a), start_angle + angle_original, normal);
						glm::vec3 direction_adjusted = glm::rotate(glm::normalize(pos_b - pos_a), start_angle + angle_adjusted, normal);

						travel_field_original_positions.push_back(origin);
						travel_field_original_positions.push_back(origin + direction_original * 0.1f);
						travel_field_original_indices.push_back(travel_field_original_positions.size() - 2);
						travel_field_original_indices.push_back(travel_field_original_positions.size() - 1);

						travel_field_adjusted_positions.push_back(origin);
						travel_field_adjusted_positions.push_back(origin + direction_adjusted * 0.1f);
						travel_field_adjusted_indices.push_back(travel_field_adjusted_positions.size() - 2);
						travel_field_adjusted_indices.push_back(travel_field_adjusted_positions.size() - 1);
					}
					mesh_buffer_travel_field_original.load(travel_field_original_positions, travel_field_original_indices);
					mesh_buffer_travel_field_adjusted.load(travel_field_adjusted_positions, travel_field_adjusted_indices);
				}
			}
			else {
				// Create field data for rendering
				{
					std::vector<glm::vec3> field_positions;
					std::vector<int> field_indices;

					for (int f_idx = 0; f_idx < mesh_ex.faces.size(); f_idx++) {
						const FaceEx& f = mesh_ex.faces[f_idx];
						glm::vec3 normal = f.normal;
						float angle = field_angles[f_idx];

						glm::vec3 origin = mesh_ex.centerOfMass(f_idx);
						glm::vec3 direction = glm::rotate(glm::normalize(mesh_ex.vertexToVertex(f.vertices[0], f.vertices[1])), start_angle + angle, normal);

						field_positions.push_back(origin);
						field_positions.push_back(origin + direction * 0.1f);
						field_indices.push_back(field_positions.size() - 2);
						field_indices.push_back(field_positions.size() - 1);
					}
					mesh_buffer_field.load(field_positions, field_indices);
				}
			}
		}

		// Reset screen buffers
		glViewport(0, 0, window.getWindowSize().x, window.getWindowSize().y);
		glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Get camera data
		const glm::vec3 cameraPos = trackball.position();
		const glm::mat4 model{ 1.0f };
		const glm::mat4 view = trackball.viewMatrix();
		const glm::mat4 projection = trackball.projectionMatrix();
		const glm::mat4 mvp = projection * view * model;

		{
			// Draw mesh
			normal_shader.bind();
			glUniformMatrix4fv(normal_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));

			glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_main.vbo);
			glBindVertexArray(mesh_buffer_main.vao);
			glVertexAttribPointer(normal_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
			glVertexAttribPointer(normal_shader.getAttributeLocation("normal"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
			glDrawElements(GL_TRIANGLES, mesh_buffer_main.indices_amount, GL_UNSIGNED_INT, nullptr);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindVertexArray(0);
		}

		std::vector<SingularityUniform> singularity_uniforms = {};
		int singularity_amount = 0;
		GLuint singularity_ubo;
		if (stage == VERTEX_K || stage >= INSPECT_RESULT) {
			// Set vertex uniform data for instanced singularity drawing
			for (int v_idx = 0; v_idx < mesh_ex.vertices.size(); v_idx++)
				singularity_uniforms.push_back(SingularityUniform::fromVertexEx(mesh_ex, v_idx));
			singularity_amount = singularity_uniforms.size();

			glGenBuffers(1, &singularity_ubo);
			glBindBuffer(GL_UNIFORM_BUFFER, singularity_ubo);
			glBufferData(GL_UNIFORM_BUFFER, singularity_amount * sizeof(SingularityUniform) + 16, NULL, GL_DYNAMIC_DRAW);
			glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(int), &singularity_amount);
			glBufferSubData(GL_UNIFORM_BUFFER, 16, singularity_amount * sizeof(SingularityUniform), singularity_uniforms.data());
			glBindBuffer(GL_UNIFORM_BUFFER, 0);
		}

		// Draw background
		{
			glDisable(GL_DEPTH_TEST);
			glDepthMask(GL_FALSE);

			// Draw wireframe background
			glm::vec4 dark_grey_transparent{ 0.3f, 0.3f, 0.3f, 0.15f };
			wireframe_shader.bind();
			glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
			glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(dark_grey_transparent));
			glLineWidth(1);

			glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_wireframe.vbo);
			glBindVertexArray(mesh_buffer_wireframe.vao);
			glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
			glDrawElements(GL_LINES, mesh_buffer_wireframe.indices_amount, GL_UNSIGNED_INT, nullptr);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindVertexArray(0);

			if (stage == INSPECT_RESULT) {
				if (dual_edges_path.size()) {
					if (show_travel_original) {
						// Draw travel vector (original)
						wireframe_shader.bind();
						glm::vec4 red_transparent{ 1.0f, 0.2f, 0.2f, 0.35f };
						glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
						glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(red_transparent));
						glLineWidth(8);

						glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_travel_vector_original.vbo);
						glBindVertexArray(mesh_buffer_travel_vector_original.vao);
						glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
						glDrawElements(GL_LINES, mesh_buffer_travel_vector_original.indices_amount, GL_UNSIGNED_INT, 0);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
						glBindVertexArray(0);
					}

					if (show_travel_adjusted) {
						// Draw travel vector (adjusted)
						wireframe_shader.bind();
						glm::vec4 blue_transparent{ 0.1f, 0.4f, 1.0f, 0.35f };
						glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
						glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(blue_transparent));
						glLineWidth(8);

						glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_travel_vector_adjusted.vbo);
						glBindVertexArray(mesh_buffer_travel_vector_adjusted.vao);
						glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
						glDrawElements(GL_LINES, mesh_buffer_travel_vector_adjusted.indices_amount, GL_UNSIGNED_INT, 0);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
						glBindVertexArray(0);
					}
				}
			}

			if (stage == VERTEX_K) {
				// Draw background singularity vertices
				singularity_shader.bind();
				glUniformMatrix4fv(singularity_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
				glUniform1i(singularity_shader.getUniformLocation("selected_vertex_idx"), selected_vertex_idx);
				glUniform1i(singularity_shader.getUniformLocation("draw_zero_k"), stage == VERTEX_K);
				glUniform1f(singularity_shader.getUniformLocation("alpha"), 0.35f);
				singularity_shader.bindUniformBlock("vertex_buffer", 0, singularity_ubo);

				glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_singularity.vbo);
				glBindVertexArray(mesh_buffer_singularity.vao);
				glVertexAttribPointer(singularity_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
				glDrawElementsInstanced(GL_TRIANGLES, mesh_buffer_singularity.indices_amount, GL_UNSIGNED_INT, 0, singularity_amount);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
				glBindVertexArray(0);
			}

			if (stage == TREE_COTREE) {
				if (show_noncons) {
					// Draw selected noncon-cycle
					int selected_noncon_idx = selected_noncon_item_idx - 1;
					if (selected_noncon_idx > -1) {
						const MeshBuffer& mesh_buffer_noncon_selected = mesh_buffers_noncon[selected_noncon_idx];

						int k = noncon_ks[selected_noncon_idx];
						float a = pow(0.6, abs(k));
						float b = pow(0.9, abs(k));
						glm::vec4 noncon_color = k > 0 ? glm::vec4{ 1.0f, b, a, 0.35f } : glm::vec4{ a, b, 1.0f, 0.35f };

						wireframe_shader.bind();
						glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
						glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(noncon_color));
						glLineWidth(20);

						glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_noncon_selected.vbo);
						glBindVertexArray(mesh_buffer_noncon_selected.vao);
						glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
						glDrawElements(GL_LINES, mesh_buffer_noncon_selected.indices_amount, GL_UNSIGNED_INT, 0);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
						glBindVertexArray(0);
					}
				}
			}

			// Draw cycle background
			if (stage == INSPECT_RESULT) {
				wireframe_shader.bind();
				glm::vec4 red_transparent{ 1.0f, 0.2f, 0.2f, 0.35f };
				glm::vec4 green_transparent{ 0.2f, 1.0f, 0.2f, 0.35f };
				glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
				glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(dual_edges_path.empty() ? red_transparent : green_transparent));
				glLineWidth(5);

				glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_dual_edges.vbo);
				glBindVertexArray(mesh_buffer_dual_edges.vao);
				glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
				glDrawElements(GL_LINES, mesh_buffer_dual_edges.indices_amount, GL_UNSIGNED_INT, 0);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
				glBindVertexArray(0);
			}

			glEnable(GL_DEPTH_TEST);
			glDepthMask(GL_TRUE);
		}

		// Draw foreground
		{
			if (stage == VERTEX_K || (stage == INSPECT_RESULT && show_singularity_data)) {
				float alpha = stage == VERTEX_K ? 1.0f : 0.5f;

				// Draw singularity vertices
				singularity_shader.bind();
				glUniformMatrix4fv(singularity_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
				glUniform1i(singularity_shader.getUniformLocation("selected_vertex_idx"), selected_vertex_idx);
				glUniform1i(singularity_shader.getUniformLocation("draw_zero_k"), stage == VERTEX_K);
				glUniform1f(singularity_shader.getUniformLocation("alpha"), alpha);
				singularity_shader.bindUniformBlock("vertex_buffer", 0, singularity_ubo);

				glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_singularity.vbo);
				glBindVertexArray(mesh_buffer_singularity.vao);
				glVertexAttribPointer(singularity_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
				glDrawElementsInstanced(GL_TRIANGLES, mesh_buffer_singularity.indices_amount, GL_UNSIGNED_INT, 0, singularity_amount);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
				glBindVertexArray(0);
			}

			if ((stage == TREE_COTREE && show_noncons) || (stage == INSPECT_RESULT && show_singularity_data)) {
				int selected_noncon_idx = -1;

				if (stage == TREE_COTREE) {
					// Draw selected noncon-cycle
					selected_noncon_idx = selected_noncon_item_idx - 1;
					if (selected_noncon_idx > -1) {
						const MeshBuffer& mesh_buffer_noncon_selected = mesh_buffers_noncon[selected_noncon_idx];

						int k = noncon_ks[selected_noncon_idx];
						float a = pow(0.6, abs(k));
						float b = pow(0.9, abs(k));
						glm::vec4 noncon_color = k > 0 ? glm::vec4{ 1.0f, b, a, 1.0f } : glm::vec4{ a, b, 1.0f, 1.0f };

						wireframe_shader.bind();
						glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
						glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(noncon_color));
						glLineWidth(20);

						glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_noncon_selected.vbo);
						glBindVertexArray(mesh_buffer_noncon_selected.vao);
						glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
						glDrawElements(GL_LINES, mesh_buffer_noncon_selected.indices_amount, GL_UNSIGNED_INT, 0);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
						glBindVertexArray(0);
					}
				}

				// Draw unselected noncon-cycles
				for (int i = 0; i < mesh_buffers_noncon.size(); i++) {
					if (i == selected_noncon_idx)
						continue;
					const MeshBuffer& mesh_buffer_noncon = mesh_buffers_noncon[i];

					int k = noncon_ks[i];
					float a = pow(0.6, abs(k));
					float b = pow(0.9, abs(k));
					if (stage == INSPECT_RESULT && k == 0)
						continue;
					float alpha = stage == TREE_COTREE ? 1.0f : 0.5f;
					glm::vec4 noncon_color = k > 0 ? glm::vec4{ 1.0f, b, a, alpha } : glm::vec4{ a, b, 1.0f, alpha };

					wireframe_shader.bind();
					glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
					glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(noncon_color));
					glLineWidth(12);

					glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_noncon.vbo);
					glBindVertexArray(mesh_buffer_noncon.vao);
					glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
					glDrawElements(GL_LINES, mesh_buffer_noncon.indices_amount, GL_UNSIGNED_INT, 0);
					glBindBuffer(GL_ARRAY_BUFFER, 0);
					glBindVertexArray(0);
				}
			}

			if (stage == INSPECT_RESULT) {
				if (!dual_edges_path.empty()) {
					if (show_travel_original) {
						// Draw travel vector (original)
						wireframe_shader.bind();
						glm::vec4 red{ 1.0f, 0.2f, 0.2f, 1.0f };
						glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
						glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(red));
						glLineWidth(8);

						glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_travel_vector_original.vbo);
						glBindVertexArray(mesh_buffer_travel_vector_original.vao);
						glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
						glDrawElements(GL_LINES, mesh_buffer_travel_vector_original.indices_amount, GL_UNSIGNED_INT, 0);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
						glBindVertexArray(0);
					}

					if (show_travel_adjusted) {
						// Draw travel vector (adjusted)
						wireframe_shader.bind();
						glm::vec4 blue{ 0.1f, 0.4f, 1.0f, 1.0f };
						glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
						glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(blue));
						glLineWidth(8);

						glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_travel_vector_adjusted.vbo);
						glBindVertexArray(mesh_buffer_travel_vector_adjusted.vao);
						glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
						glDrawElements(GL_LINES, mesh_buffer_travel_vector_adjusted.indices_amount, GL_UNSIGNED_INT, 0);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
						glBindVertexArray(0);
					}

					if (show_travel_original) {
						// Draw travel field (original)
						wireframe_shader.bind();
						glm::vec4 red_transparent{ 1.0f, 0.2f, 0.2f, 0.35f };
						glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
						glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(red_transparent));
						glLineWidth(5);

						glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_travel_field_original.vbo);
						glBindVertexArray(mesh_buffer_travel_field_original.vao);
						glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
						glDrawElements(GL_LINES, mesh_buffer_travel_field_original.indices_amount, GL_UNSIGNED_INT, 0);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
						glBindVertexArray(0);
					}

					if (show_travel_adjusted) {
						// Draw travel field (adjusted)
						wireframe_shader.bind();
						glm::vec4 blue_transparent{ 0.1f, 0.4f, 1.0f, 0.35f };
						glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
						glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(blue_transparent));
						glLineWidth(5);

						glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_travel_field_adjusted.vbo);
						glBindVertexArray(mesh_buffer_travel_field_adjusted.vao);
						glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
						glDrawElements(GL_LINES, mesh_buffer_travel_field_adjusted.indices_amount, GL_UNSIGNED_INT, 0);
						glBindBuffer(GL_ARRAY_BUFFER, 0);
						glBindVertexArray(0);
					}
				}
				else {
					// Draw field
					wireframe_shader.bind();
					glm::vec4 blue_transparent{ 0.1f, 0.4f, 1.0f, 1.0f };
					glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
					glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(blue_transparent));
					glLineWidth(5);
					glPointSize(10);

					glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_field.vbo);
					glBindVertexArray(mesh_buffer_field.vao);
					glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
					glDrawElements(GL_LINES, mesh_buffer_field.indices_amount, GL_UNSIGNED_INT, 0);
					glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3) * 2, (void*)0);
					glDrawElements(GL_POINTS, mesh_buffer_field.indices_amount / 2, GL_UNSIGNED_INT, 0);
					glBindBuffer(GL_ARRAY_BUFFER, 0);
					glBindVertexArray(0);
				}
			}

			if (stage == TREE_COTREE) {
				// Draw tree
				if (show_tree) {
					wireframe_shader.bind();
					glm::vec4 green{ 0.2f, 1.0f, 0.2f, 1.0f };
					glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
					glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(green));
					glLineWidth(8);

					glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_tree.vbo);
					glBindVertexArray(mesh_buffer_tree.vao);
					glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
					glDrawElements(GL_LINES, mesh_buffer_tree.indices_amount, GL_UNSIGNED_INT, 0);
					glBindBuffer(GL_ARRAY_BUFFER, 0);
					glBindVertexArray(0);
				}

				if (show_cotree) {
					// Draw cotree
					wireframe_shader.bind();
					glm::vec4 magenta{ 0.8f, 0.2f, 1.0f, 1.0f };
					glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
					glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(magenta));
					glLineWidth(8);

					glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_cotree.vbo);
					glBindVertexArray(mesh_buffer_cotree.vao);
					glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
					glDrawElements(GL_LINES, mesh_buffer_cotree.indices_amount, GL_UNSIGNED_INT, 0);
					glBindBuffer(GL_ARRAY_BUFFER, 0);
					glBindVertexArray(0);
				}
			}

			if (stage == INSPECT_RESULT) {
				// Draw cycle
				wireframe_shader.bind();
				glm::vec4 red{ 1.0f, 0.2f, 0.2f, 1.0f };
				glm::vec4 green{ 0.2f, 1.0f, 0.2f, 1.0f };
				glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
				glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(dual_edges_path.empty() ? red : green));
				glLineWidth(12);

				glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_dual_edges.vbo);
				glBindVertexArray(mesh_buffer_dual_edges.vao);
				glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
				glDrawElements(GL_LINES, mesh_buffer_dual_edges.indices_amount, GL_UNSIGNED_INT, 0);
				glBindBuffer(GL_ARRAY_BUFFER, 0);
				glBindVertexArray(0);
			}

			// Draw wireframe
			glm::vec4 dark_grey{ 0.3f, 0.3f, 0.3f, 1.0f };
			wireframe_shader.bind();
			glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
			glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(dark_grey));
			glLineWidth(2);

			glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_wireframe.vbo);
			glBindVertexArray(mesh_buffer_wireframe.vao);
			glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
			glDrawElements(GL_LINES, mesh_buffer_wireframe.indices_amount, GL_UNSIGNED_INT, nullptr);
			glBindBuffer(GL_ARRAY_BUFFER, 0);
			glBindVertexArray(0);
		}

		window.swapBuffers();
	}

	mesh_buffer_main.cleanUp();
	mesh_buffer_wireframe.cleanUp();
	mesh_buffer_singularity.cleanUp();

	mesh_buffer_tree.cleanUp();
	mesh_buffer_cotree.cleanUp();
	for (const MeshBuffer& mb : mesh_buffers_noncon)
		mb.cleanUp();

	return 0;
}
