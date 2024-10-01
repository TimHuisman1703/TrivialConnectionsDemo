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
#include <vector>
#include <array>

#include "buffer.h"
#include "mesh_ex.h"

// Configuration
const int WIDTH = 1200;
const int HEIGHT = 800;

struct Texture {
    int width;
    int height;
    int channels;
    stbi_uc* texture_data;
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
    GLuint vbo_dual_edges;
    glGenBuffers(1, &vbo_dual_edges);
    GLuint ibo_dual_edges;
    glGenBuffers(1, &ibo_dual_edges);
    GLuint vao_dual_edges;
    glGenVertexArrays(1, &vao_dual_edges);

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
    MeshBuffer mesh_buffer = std::get<0>(mesh_data);
    MeshBuffer mesh_buffer_wireframe = std::get<1>(mesh_data);
    MeshEx mesh_ex = std::get<2>(mesh_data);

    // Load singularity mesh
    const Mesh mesh_singularity = loadMesh(std::string(RESOURCE_ROOT) + "resources/singularity.obj")[0];
    const MeshBuffer mesh_buffer_singularity = MeshBuffer::fromTriangles(mesh_singularity.vertices, mesh_singularity.triangles);

    // Parameters
    bool show_imgui = true;
    int stage = 0;
    int selected_vertex_idx = 0;
    std::vector<bool> dual_edge_selected;

    // Render data
    bool reload_dual_edges = true;
    std::vector<glm::vec3> positions_dual_edges;
    std::vector<int> indices_dual_edges;
    std::vector<int> path_dual_edges;

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

        if (stage == 1) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
                const auto opt_world_point = getWorldPositionOfPixel(trackball, window.getCursorPixel());
                if (opt_world_point)
                    selected_vertex_idx = getClosestVertexIndex(mesh_ex, *opt_world_point);
                return;
            }
        } else if (stage == 2) {
            if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
                const auto opt_world_point = getWorldPositionOfPixel(trackball, window.getCursorPixel());
                if (opt_world_point) {
                    int closest_edge_index = getClosestEdgeIndex(mesh_ex, *opt_world_point);
                    dual_edge_selected[closest_edge_index] = !dual_edge_selected[closest_edge_index];
                    reload_dual_edges = true;
                }
                return;
            }
        }
        });

    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_BLEND);

    const std::array<std::string, 3> stage_names{
        "Choosing Mesh",
        "Setting Singularities",
        "Test Cycles",
    };
    const std::array<std::string, 3> stage_descriptions{
        "Select a mesh to construct a\n  vector field on.",
        "Set the singularities of each\n  vertex.\nSelect a vertex by holding\n  SHIFT and clicking on the\n  mesh.\nIncrease or decrease its index\n  in the GUI.",
        "Test Cycles",
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
        bool prev_enabled = true;
        bool next_enabled = true;
        if (stage == 0) {
            prev_enabled = false;
        }
        else if (stage == 1) {
            if (total_k != euler_mesh_characteristic)
                next_enabled = false;
        }
        else if (stage == 2) {
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
                    stage--;

                    selected_vertex_idx = -1;
                    reload_dual_edges = true;
                }
            }
            ImGui::SameLine();
            if (ImGui::Button(next_enabled ? "Next" : "##", ImVec2(halfWidth, 19))) {
                if (next_enabled) {
                    stage++;

                    selected_vertex_idx = -1;
                    reload_dual_edges = true;

                    if (stage == 2) {
                        dual_edge_selected = {};
                        for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++)
                            dual_edge_selected.push_back(false);
                    }
                }
            }

            ImGui::Separator();
            ImGui::Text(stage_descriptions[stage].c_str());
            ImGui::Separator();

            if (stage == 0) {
                if (ImGui::Button("Icosphere (g = 0)", ImVec2(width, 19))) {
                    std::tuple<MeshBuffer, MeshBuffer, MeshEx> mesh_data = loadMainMesh("genus_0");
                    mesh_buffer = std::get<0>(mesh_data);
                    mesh_buffer_wireframe = std::get<1>(mesh_data);
                    mesh_ex = std::get<2>(mesh_data);
                }

                if (ImGui::Button("Thorus (g = 1)", ImVec2(width, 19))) {
                    std::tuple<MeshBuffer, MeshBuffer, MeshEx> mesh_data = loadMainMesh("genus_1");
                    mesh_buffer = std::get<0>(mesh_data);
                    mesh_buffer_wireframe = std::get<1>(mesh_data);
                    mesh_ex = std::get<2>(mesh_data);
                }

                if (ImGui::Button("Bunny (g = 0)", ImVec2(width, 19))) {
                    std::tuple<MeshBuffer, MeshBuffer, MeshEx> mesh_data = loadMainMesh("bunny");
                    mesh_buffer = std::get<0>(mesh_data);
                    mesh_buffer_wireframe = std::get<1>(mesh_data);
                    mesh_ex = std::get<2>(mesh_data);
                }
            } else if (stage == 1) {
                if (selected_vertex_idx > -1) {
                    ImGui::Text("Vertex Selected: #%i", selected_vertex_idx);
                    ImGui::InputInt("k", &mesh_ex.vertices[selected_vertex_idx].k);
                }
                else {
                    ImGui::Text("Vertex Selected: None");
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
            else if (stage == 2) {
                if (path_dual_edges.size() > 0) {
                    ImGui::TextColored(ImVec4(0.2, 1.0, 0.4, 1.0), "Cycle complete");
                    ImGui::Text("  #edges = %i", path_dual_edges.size());
                }
                else {
                    ImGui::TextColored(ImVec4(1.0, 0.4, 0.2, 1.0), "Cycle incomplete");
                }
            }

            ImGui::End();
            ImGui::Render();
        }

        // Logic
        if (stage == 2) {
            if (reload_dual_edges) {
                // Trace cycle
                path_dual_edges = {};
                std::vector<std::vector<int>> edges_per_face(mesh_ex.faces.size());

                int f_start_idx = -1;
                int e_start_idx = -1;
                int expected_length = 0;
                for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++) {
                    if (dual_edge_selected[e_idx]) {
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
                    path_dual_edges = { e_start_idx };

                    int f_curr_idx = mesh_ex.otherFace(e_start_idx, f_start_idx);
                    int e_curr_idx = e_start_idx;
                    while (f_curr_idx != f_start_idx) {
                        std::vector<int> edges = edges_per_face[f_curr_idx];

                        if (edges.size() != 2) {
                            path_dual_edges = {};
                            break;
                        }

                        if (edges[0] != e_curr_idx)
                            e_curr_idx = edges[0];
                        else
                            e_curr_idx = edges[1];

                        f_curr_idx = mesh_ex.otherFace(e_curr_idx, f_curr_idx);
                        path_dual_edges.push_back(e_curr_idx);
                    }

                    if (path_dual_edges.size() != expected_length)
                        path_dual_edges = {};
                }
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

        // Draw mesh
        normal_shader.bind();
        glUniformMatrix4fv(normal_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));

        glBindBuffer(GL_ARRAY_BUFFER, mesh_buffer.vbo);
        glBindVertexArray(mesh_buffer.vao);
        glVertexAttribPointer(normal_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, position));
        glVertexAttribPointer(normal_shader.getAttributeLocation("normal"), 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)offsetof(Vertex, normal));
        glDrawElements(GL_TRIANGLES, mesh_buffer.indices_amount, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        // Draw wireframe background
        glDisable(GL_DEPTH_TEST);
        glDepthMask(0.0f);

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

        glEnable(GL_DEPTH_TEST);

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

        glDepthMask(1.0f);

        // Set vertex uniform data for instanced singularity drawing
        std::vector<SingularityUniform> singularity_uniforms = {};
        for (int v_idx = 0; v_idx < mesh_ex.vertices.size(); v_idx++)
            singularity_uniforms.push_back(SingularityUniform::fromVertexEx(mesh_ex, v_idx));
        int singularity_amount = singularity_uniforms.size();

        GLuint singularity_ubo;
        glGenBuffers(1, &singularity_ubo);
        glBindBuffer(GL_UNIFORM_BUFFER, singularity_ubo);
        glBufferData(GL_UNIFORM_BUFFER, singularity_amount * sizeof(SingularityUniform) + 16, NULL, GL_DYNAMIC_DRAW);
        glBufferSubData(GL_UNIFORM_BUFFER, 0, sizeof(int), &singularity_amount);
        glBufferSubData(GL_UNIFORM_BUFFER, 16, singularity_amount * sizeof(SingularityUniform), singularity_uniforms.data());
        glBindBuffer(GL_UNIFORM_BUFFER, 0);

        // Draw background singularity vertices
        glDisable(GL_DEPTH_TEST);

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

        glEnable(GL_DEPTH_TEST);

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

        if (stage == 2) {
            if (reload_dual_edges) {
                // Create cycle data
                positions_dual_edges = {};
                indices_dual_edges = {};

                for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++) {
                    if (dual_edge_selected[e_idx]) {
                        const EdgeEx& e = mesh_ex.edges[e_idx];

                        glm::vec3 e_center = 0.5f * (mesh_ex.vertices[e.vertices[0]].position + mesh_ex.vertices[e.vertices[1]].position);
                        glm::vec3 f_a_center = mesh_ex.centerOfMass(e.faces[0]);
                        glm::vec3 f_b_center = mesh_ex.centerOfMass(e.faces[1]);

                        positions_dual_edges.push_back(f_a_center);
                        positions_dual_edges.push_back(e_center);
                        positions_dual_edges.push_back(f_b_center);

                        indices_dual_edges.push_back(positions_dual_edges.size() - 3);
                        indices_dual_edges.push_back(positions_dual_edges.size() - 2);
                        indices_dual_edges.push_back(positions_dual_edges.size() - 2);
                        indices_dual_edges.push_back(positions_dual_edges.size() - 1);
                    }
                }

                glBindBuffer(GL_ARRAY_BUFFER, vbo_dual_edges);
                glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(positions_dual_edges.size() * sizeof(glm::vec3)), positions_dual_edges.data(), GL_STATIC_DRAW);
                glBindBuffer(GL_ARRAY_BUFFER, 0);

                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_dual_edges);
                glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(indices_dual_edges.size() * sizeof(decltype(indices_dual_edges)::value_type)), indices_dual_edges.data(), GL_STATIC_DRAW);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

                glBindVertexArray(vao_dual_edges);
                glBindBuffer(GL_ARRAY_BUFFER, vbo_dual_edges);
                glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo_dual_edges);
                glEnableVertexAttribArray(0);
                glBindVertexArray(0);

                reload_dual_edges = false;
            }

            // Draw cycle background
            glDisable(GL_DEPTH_TEST);

            wireframe_shader.bind();
            glm::vec4 red_transparent{ 1.0f, 0.2f, 0.2f, 0.35f };
            glm::vec4 green_transparent{ 0.2f, 1.0f, 0.2f, 0.35f };
            glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
            glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(path_dual_edges.empty() ? red_transparent : green_transparent));
            glLineWidth(5);

            glBindBuffer(GL_ARRAY_BUFFER, vbo_dual_edges);
            glBindVertexArray(vao_dual_edges);
            glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
            glDrawElements(GL_LINES, static_cast<GLsizei>(3 * indices_dual_edges.size()), GL_UNSIGNED_INT, 0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindVertexArray(0);

            glEnable(GL_DEPTH_TEST);

            // Draw cycle
            wireframe_shader.bind();
            glm::vec4 red{ 1.0f, 0.2f, 0.2f, 1.0f };
            glm::vec4 green{ 0.2f, 1.0f, 0.2f, 1.0f };
            glUniformMatrix4fv(wireframe_shader.getUniformLocation("mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
            glUniform4fv(wireframe_shader.getUniformLocation("albedo"), 1, glm::value_ptr(path_dual_edges.empty() ? red : green));
            glLineWidth(12);

            glBindBuffer(GL_ARRAY_BUFFER, vbo_dual_edges);
            glBindVertexArray(vao_dual_edges);
            glVertexAttribPointer(wireframe_shader.getAttributeLocation("pos"), 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), (void*)0);
            glDrawElements(GL_LINES, static_cast<GLsizei>(3 * indices_dual_edges.size()), GL_UNSIGNED_INT, 0);
            glBindBuffer(GL_ARRAY_BUFFER, 0);
            glBindVertexArray(0);
        }

        window.swapBuffers();
    }

    mesh_buffer.cleanUp();
    mesh_buffer_wireframe.cleanUp();
    mesh_buffer_singularity.cleanUp();

    return 0;
}
