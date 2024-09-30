// Disable compiler warnings in third-party code (which we cannot change).
#include <framework/disable_all_warnings.h>
#include <framework/opengl_includes.h>
DISABLE_WARNINGS_PUSH()
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
#include <cstdlib> // EXIT_FAILURE
#include <framework/mesh.h>
#include <framework/shader.h>
#include <framework/trackball.h>
#include <framework/window.h>
#include <imgui/imgui.h>
#include <imgui/imgui_impl_glfw.h>
#include <imgui/imgui_impl_opengl2.h>

#include "mesh_ex.h";

struct MeshBuffer {
	GLuint vbo;
	GLuint ibo;
	GLuint vao;

	int indices_amount;
	std::vector<Vertex> vertices;
	std::vector<glm::uvec2> edges;
	std::vector<glm::uvec3> triangles;

	static MeshBuffer fromEdges(std::vector<Vertex> vertices, std::vector<glm::uvec2> edges);
	static MeshBuffer fromTriangles(std::vector<Vertex> vertices, std::vector<glm::uvec3> triangles);

	void cleanUp() const;
};

struct SingularityUniform {
	glm::vec3 position;
	int k;

	static SingularityUniform fromVertexEx(const MeshEx& mesh_ex, int v_idx);
};
