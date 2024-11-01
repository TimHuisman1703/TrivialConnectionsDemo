#include "buffer.h"

MeshBuffer MeshBuffer::empty()
{
	MeshBuffer result;

	glGenBuffers(1, &result.vbo);
	glGenBuffers(1, &result.ibo);
	glGenVertexArrays(1, &result.vao);

	return result;
}

MeshBuffer MeshBuffer::fromEdges(std::vector<Vertex> vertices, std::vector<glm::uvec2> edges)
{
	MeshBuffer result;

	glGenBuffers(1, &result.vbo);
	glBindBuffer(GL_ARRAY_BUFFER, result.vbo);
	glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(vertices.size() * sizeof(Vertex)), vertices.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glGenBuffers(1, &result.ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, result.ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(edges.size() * sizeof(decltype(edges)::value_type)), edges.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	glGenVertexArrays(1, &result.vao);
	glBindVertexArray(result.vao);
	glBindBuffer(GL_ARRAY_BUFFER, result.vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, result.ibo);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glBindVertexArray(0);

	result.indices_amount = static_cast<GLsizei>(2 * edges.size());
	result.vertices = vertices;
	result.edges = edges;

	return result;
}

MeshBuffer MeshBuffer::fromTriangles(std::vector<Vertex> vertices, std::vector<glm::uvec3> triangles)
{
	MeshBuffer result;

	glGenBuffers(1, &result.vbo);
	glBindBuffer(GL_ARRAY_BUFFER, result.vbo);
	glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(vertices.size() * sizeof(Vertex)), vertices.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glGenBuffers(1, &result.ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, result.ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(triangles.size() * sizeof(decltype(triangles)::value_type)), triangles.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	glGenVertexArrays(1, &result.vao);
	glBindVertexArray(result.vao);
	glBindBuffer(GL_ARRAY_BUFFER, result.vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, result.ibo);
	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glBindVertexArray(0);

	result.indices_amount = static_cast<GLsizei>(3 * triangles.size());
	result.vertices = vertices;
	result.triangles = triangles;

	return result;
}

void MeshBuffer::load(std::vector<glm::vec3> positions, std::vector<int> indices) {
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(positions.size() * sizeof(glm::vec3)), positions.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, static_cast<GLsizeiptr>(indices.size() * sizeof(decltype(indices)::value_type)), indices.data(), GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	glBindVertexArray(vao);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
	glEnableVertexAttribArray(0);
	glBindVertexArray(0);

	indices_amount = static_cast<GLsizei>(indices.size());
}

void MeshBuffer::cleanUp() const {
	glDeleteBuffers(1, &vbo);
	glDeleteBuffers(1, &ibo);
	glDeleteVertexArrays(1, &vao);
}

SingularityUniform SingularityUniform::fromVertexEx(const MeshEx& mesh_ex, int v_idx) {
	const VertexEx& v = mesh_ex.vertices[v_idx];
	return { v.position, v.k };
}
