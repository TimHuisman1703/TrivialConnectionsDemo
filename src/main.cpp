// Disable compiler warnings in third-party code (which we cannot change).
#include <framework/disable_all_warnings.h>
#include <framework/opengl_includes.h>
DISABLE_WARNINGS_PUSH()
#include <fmt/format.h>
// Include glad before glfw3
#include <GLFW/glfw3.h>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
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
#include "tree_cotree.h"

// Configuration
const int WIDTH = 1200;
const int HEIGHT = 800;

struct Texture {
    int width;
    int height;
    int channels;
    stbi_uc* texture_data;
};

enum Stage {
    CHOOSE_MESH,
    VERTEX_K,
    TREE_COTREE,
    TEST_CYCLES,

    STAGE_NR_ITEMS
};

std::tuple<MeshBuffer, MeshBuffer, MeshEx> loadMainMesh(std::string name) {
    // Load main mesh
    const Mesh mesh = loadMesh(std::string(RESOURCE_ROOT) + fmt::format("resources/{}.obj", name.c_str()))[0];
    const MeshBuffer mesh_buffer = MeshBuffer::fromTriangles(mesh.vertices, mesh.triangles);

    // Create extended mesh representation
    MeshEx mesh_ex = MeshEx::fromMesh(mesh);
    std::vector<glm::uvec2> mesh_edges = {};
    for (const EdgeEx& e : mesh_ex.edges) {
        mesh_edges.push_back(e.vertices);
    }

    // Create wireframe buffer
    const MeshBuffer mesh_buffer_wireframe = MeshBuffer::fromEdges(mesh.vertices, mesh_edges);

    return { mesh_buffer, mesh_buffer_wireframe, mesh_ex };
}

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
    Window window { "Trivial Connections - Demonstration", glm::ivec2(WIDTH, HEIGHT), OpenGLVersion::GL45 };

    // Camera settings
    glm::vec3 look_at{ 0.0f };
    glm::vec3 rotations{ 0.6f, 0.3f, 0.0f };
    float fov_y = 30.0f;
    float dist = 4.5f;
    Trackball trackball { &window, glm::radians(fov_y) };
    trackball.setCamera(look_at, rotations, dist);

    // Create buffers
    MeshBuffer mesh_buffer_dual_edges = MeshBuffer::empty();
    MeshBuffer mesh_buffer_tree = MeshBuffer::empty();
    MeshBuffer mesh_buffer_cotree = MeshBuffer::empty();
    MeshBuffer mesh_buffer_noncon = MeshBuffer::empty();

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

    std::tuple<MeshBuffer, MeshBuffer, MeshEx> mesh_data = loadMainMesh("genus_1");
    MeshBuffer mesh_buffer_object = std::get<0>(mesh_data);
    MeshBuffer mesh_buffer_wireframe = std::get<1>(mesh_data);
    MeshEx mesh_ex = std::get<2>(mesh_data);

    // Load singularity mesh
    const Mesh mesh_singularity = loadMesh(std::string(RESOURCE_ROOT) + "resources/singularity.obj")[0];
    const MeshBuffer mesh_buffer_singularity = MeshBuffer::fromTriangles(mesh_singularity.vertices, mesh_singularity.triangles);

    // Parameters
    bool show_imgui = true;
    Stage stage = CHOOSE_MESH;
    int selected_vertex_idx = 0;

    // Tree-cotree data
    std::vector<int> tree_assignment;

    // Cycle data
    std::vector<bool> dual_edges_selected;
    bool dual_edges_reload = true;
    std::vector<int> dual_edges_path;
    std::vector<bool> dual_edges_vertex_partition;
    double curvature_partition = 0.0;

    // Handle key press
    window.registerKeyCallback([&](int key, int scancode, int action, int mods) {
        if (key == '\\' && action == GLFW_PRESS) {
            show_imgui = !show_imgui;
        }

        if (action != GLFW_RELEASE)
            return;
        });

    // Handle mouse 
    window.registerMouseButtonCallback([&](int button, int action, int mods) {
        if (!window.isKeyPressed(GLFW_KEY_LEFT_SHIFT))
            return;

        if (stage == VERTEX_K) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
                const auto opt_world_point = getWorldPositionOfPixel(trackball, window.getCursorPixel());
                if (opt_world_point)
                    selected_vertex_idx = getClosestVertexIndex(mesh_ex, *opt_world_point);
                return;
            }
        } else if (stage == TEST_CYCLES) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
                const auto opt_world_point = getWorldPositionOfPixel(trackball, window.getCursorPixel());
                if (opt_world_point) {
                    int closest_edge_index = getClosestEdgeIndex(mesh_ex, *opt_world_point);
                    dual_edges_selected[closest_edge_index] = !dual_edges_selected[closest_edge_index];
                    dual_edges_reload = true;
                }
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
        "Test Cycles",
    };
    const std::array<std::string, 4> stage_descriptions{
        "Select a mesh to construct a\n  vector field on.",
        "Set the singularities of each\n  vertex.\nSelect a vertex by holding\n  SHIFT and clicking on the\n  mesh.\nIncrease or decrease its index\n  in the GUI.",
        "Tree-Cotree Decomposition (TODO)",
        "Test Cycles (TODO)",
    };

    while (!window.shouldClose()) {
        // Update input
        window.updateInput();

        // Poincare-Hopf calculations
        int total_k = 0;
        for (const VertexEx& v : mesh_ex.vertices)
            total_k += v.k;
        int euler_mesh_characteristic = mesh_ex.vertices.size() - mesh_ex.edges.size() + mesh_ex.faces.size();

        // Stage checks
        bool prev_enabled = stage > 0;
        bool next_enabled = stage < STAGE_NR_ITEMS - 1;
        if (stage == VERTEX_K) {
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
            if (ImGui::Button(next_enabled ? "Next" : "##", ImVec2(halfWidth, 19))) {
                if (next_enabled) {
                    stage = static_cast<Stage>(stage + 1);

                    selected_vertex_idx = -1;
                    dual_edges_reload = true;

                    if (stage == TREE_COTREE) {
                        tree_assignment = treeCotreeDecompose(mesh_ex);
                        std::vector<std::vector<int>> paths = findNoncontractibleCycles(mesh_ex, tree_assignment);

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
                        std::vector<glm::vec3> noncon_positions;
                        std::vector<int> noncon_indices;
                        for (std::vector<int>& path : paths) {
                            for (int e_idx : path) {
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
                        }
                        mesh_buffer_noncon.load(noncon_positions, noncon_indices);
                    }

                    if (stage == TEST_CYCLES) {
                        dual_edges_selected = {};
                        for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++)
                            dual_edges_selected.push_back(false);
                        dual_edges_vertex_partition = {};
                        for (int v_idx = 0; v_idx < mesh_ex.vertices.size(); v_idx++)
                            dual_edges_vertex_partition.push_back(false);
                    }
                }
            }

            ImGui::Separator();
            ImGui::Text(stage_descriptions[stage].c_str());
            ImGui::Separator();

            if (stage == CHOOSE_MESH) {
                if (ImGui::Button("Icosphere (g = 0)", ImVec2(width, 19))) {
                    std::tuple<MeshBuffer, MeshBuffer, MeshEx> mesh_data = loadMainMesh("genus_0");
                    mesh_buffer_object = std::get<0>(mesh_data);
                    mesh_buffer_wireframe = std::get<1>(mesh_data);
                    mesh_ex = std::get<2>(mesh_data);
                }

                if (ImGui::Button("Thorus (g = 1)", ImVec2(width, 19))) {
                    std::tuple<MeshBuffer, MeshBuffer, MeshEx> mesh_data = loadMainMesh("genus_1");
                    mesh_buffer_object = std::get<0>(mesh_data);
                    mesh_buffer_wireframe = std::get<1>(mesh_data);
                    mesh_ex = std::get<2>(mesh_data);
                }

                if (ImGui::Button("Bunny (g = 0)", ImVec2(width, 19))) {
                    std::tuple<MeshBuffer, MeshBuffer, MeshEx> mesh_data = loadMainMesh("bunny");
                    mesh_buffer_object = std::get<0>(mesh_data);
                    mesh_buffer_wireframe = std::get<1>(mesh_data);
                    mesh_ex = std::get<2>(mesh_data);
                }
            } else if (stage == VERTEX_K) {
                if (selected_vertex_idx > -1) {
                    ImGui::Text("Vertex Selected: #%i", selected_vertex_idx);
                    ImGui::InputInt("k", &mesh_ex.vertices[selected_vertex_idx].k);
                }
                else {
                    ImGui::Text("Vertex Selected: None");
                    ImGui::NewLine();
                }

                ImGui::Separator();
                if (total_k == euler_mesh_characteristic) {
                    ImGui::TextColored(ImVec4(0.2, 1.0, 0.4, 1.0), "Poincare-Hopf satisfied");
                    ImGui::Text("  sum(k) = %i == %i", total_k, euler_mesh_characteristic);
                }
                else {
                    ImGui::TextColored(ImVec4(1, 0.4, 0.2, 1.0), "Poincare-Hopf not satisfied");
                    ImGui::Text("  sum(k) = %i != % i", total_k, euler_mesh_characteristic);
                }
            }
            else if (stage == TEST_CYCLES) {
                if (!dual_edges_path.empty()) {
                    ImGui::TextColored(ImVec4(0.2, 1.0, 0.4, 1.0), "Cycle complete");
                    ImGui::Text("  #edges = %i", dual_edges_path.size());
                    ImGui::Text("  defect = %f", curvature_partition);
                }
                else {
                    ImGui::TextColored(ImVec4(1.0, 0.4, 0.2, 1.0), "Cycle incomplete");
                }
            }

            ImGui::End();
            ImGui::Render();
        }

        // Logic
        if (stage == TEST_CYCLES) {
            if (dual_edges_reload) {
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

                    std::fill(dual_edges_vertex_partition.begin(), dual_edges_vertex_partition.end(), false);
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

                        if (partition_size == dual_edges_vertex_partition.size()) {
                            // The cycle is non-contractible
                        }

                        if (2 * partition_size > dual_edges_vertex_partition.size()) {
                            for (int v_idx = 0; v_idx < dual_edges_vertex_partition.size(); v_idx++)
                                dual_edges_vertex_partition[v_idx] = !dual_edges_vertex_partition[v_idx];
                        }

                        curvature_partition = 0.0;
                        for (int v_idx = 0; v_idx < mesh_ex.vertices.size(); v_idx++) {
                            if (dual_edges_vertex_partition[v_idx])
                                curvature_partition += mesh_ex.defectAroundVertex(v_idx);
                        }
                    }
                }

                // Create cycle data for rendering
                std::vector<glm::vec3> dual_edges_positions;
                std::vector<int> dual_edges_indices;
                for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++) {
                    if (dual_edges_selected[e_idx]) {
                        const EdgeEx& e = mesh_ex.edges[e_idx];

                        glm::vec3 e_center = 0.5f * (mesh_ex.vertices[e.vertices[0]].position + mesh_ex.vertices[e.vertices[1]].position);
                        glm::vec3 f_a_center = mesh_ex.centerOfMass(e.faces[0]);
                        glm::vec3 f_b_center = mesh_ex.centerOfMass(e.faces[1]);

                        dual_edges_positions.push_back(f_a_center);
                        dual_edges_positions.push_back(e_center);
                        dual_edges_positions.push_back(f_b_center);

                        dual_edges_indices.push_back(dual_edges_positions.size() - 3);
                        dual_edges_indices.push_back(dual_edges_positions.size() - 2);
                        dual_edges_indices.push_back(dual_edges_positions.size() - 2);
                        dual_edges_indices.push_back(dual_edges_positions.size() - 1);
                    }
                }

                mesh_buffer_dual_edges.load(dual_edges_positions, dual_edges_indices);

                dual_edges_reload = false;
            }
        }

        // Reset screen buffers
        glViewport(0, 0, window.getWindowSize().x, window.getWindowSize().y);
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Get camera data
        const glm::vec3 cameraPos = trackball.position();
        const glm::mat4 model { 1.0f };
        const glm::mat4 view = trackball.viewMatrix();
        const glm::mat4 projection = trackball.projectionMatrix();
        const glm::mat4 mvp = projection * view * model;

        {
            // Draw mesh
            normal_shader.bind();
            glUniformMatrix4fv(normal_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));

            glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_object.vbo);
            glBindVertexArray(mesh_buffer_object.vao);
            glVertexAttribPointer(normal_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
            glVertexAttribPointer(normal_shader.getAttributeLocation("normal"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
            glDrawElements(GL_TRIANGLES, mesh_buffer_object.indices_amount, GL_UNSIGNED_INT, nullptr);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindVertexArray(0);
        }

        std::vector<SingularityUniform> singularity_uniforms = {};
        int singularity_amount = 0;
        GLuint singularity_ubo;
        if (stage >= VERTEX_K) {
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
            glDepthMask(0.0f);

            // Draw wireframe background
            glm::vec4 dark_grey_transparent{ 0.3f, 0.3f, 0.3f, 0.2f };
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

            if (stage >= VERTEX_K) {
                // Draw background singularity vertices
                singularity_shader.bind();
                glUniformMatrix4fv(singularity_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
                glUniform1i(singularity_shader.getUniformLocation("selected_vertex_idx"), selected_vertex_idx);
                glUniform1i(singularity_shader.getUniformLocation("stage"), stage);
                glUniform1i(singularity_shader.getUniformLocation("transparent"), true);
                singularity_shader.bindUniformBlock("vertex_buffer", 0, singularity_ubo);

                glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_singularity.vbo);
                glBindVertexArray(mesh_buffer_singularity.vao);
                glVertexAttribPointer(singularity_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
                glDrawElementsInstanced(GL_TRIANGLES, mesh_buffer_singularity.indices_amount, GL_UNSIGNED_INT, 0, singularity_amount);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glBindVertexArray(0);
            }

            // Draw cycle background
            if (stage == TEST_CYCLES) {
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
            glDepthMask(1.0f);
        }

        // Draw foreground
        {
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

            if (stage >= VERTEX_K) {
                // Draw singularity vertices
                singularity_shader.bind();
                glUniformMatrix4fv(singularity_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
                glUniform1i(singularity_shader.getUniformLocation("selected_vertex_idx"), selected_vertex_idx);
                glUniform1i(singularity_shader.getUniformLocation("stage"), stage);
                glUniform1i(singularity_shader.getUniformLocation("transparent"), false);
                singularity_shader.bindUniformBlock("vertex_buffer", 0, singularity_ubo);

                glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_singularity.vbo);
                glBindVertexArray(mesh_buffer_singularity.vao);
                glVertexAttribPointer(singularity_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
                glDrawElementsInstanced(GL_TRIANGLES, mesh_buffer_singularity.indices_amount, GL_UNSIGNED_INT, 0, singularity_amount);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glBindVertexArray(0);
            }

            if (stage == TREE_COTREE) {
                // Draw noncon-cycles
                wireframe_shader.bind();
                glm::vec4 green{ 0.2f, 1.0f, 0.2f, 1.0f };
                glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
                glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(green));
                glLineWidth(12);

                glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_noncon.vbo);
                glBindVertexArray(mesh_buffer_noncon.vao);
                glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
                glDrawElements(GL_LINES, mesh_buffer_noncon.indices_amount, GL_UNSIGNED_INT, 0);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glBindVertexArray(0);

                // Draw tree
                wireframe_shader.bind();
                glm::vec4 orange{ 1.0f, 0.4f, 0.2f, 1.0f };
                glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
                glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(orange));
                glLineWidth(8);

                glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_tree.vbo);
                glBindVertexArray(mesh_buffer_tree.vao);
                glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
                glDrawElements(GL_LINES, mesh_buffer_tree.indices_amount, GL_UNSIGNED_INT, 0);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glBindVertexArray(0);

                // Draw cotree
                wireframe_shader.bind();
                glm::vec4 blue{ 0.2f, 0.2f, 1.0f, 1.0f };
                glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
                glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(blue));
                glLineWidth(8);

                glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer_cotree.vbo);
                glBindVertexArray(mesh_buffer_cotree.vao);
                glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
                glDrawElements(GL_LINES, mesh_buffer_cotree.indices_amount, GL_UNSIGNED_INT, 0);
                glBindBuffer(GL_ARRAY_BUFFER, 0);
                glBindVertexArray(0);
            }

            if (stage == TEST_CYCLES) {
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
        }

        window.swapBuffers();
    }

    mesh_buffer_object.cleanUp();
    mesh_buffer_wireframe.cleanUp();
    mesh_buffer_singularity.cleanUp();

    mesh_buffer_tree.cleanUp();
    mesh_buffer_cotree.cleanUp();
    mesh_buffer_noncon.cleanUp();

    return 0;
}
