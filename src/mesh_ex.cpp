#include "mesh_ex.h"

MeshEx MeshEx::fromMesh(const Mesh mesh)
{
	MeshEx result;

	// initialize vertex positions + normals
	for (const Vertex& v_old : mesh.vertices) {
		result.vertices.push_back({ v_old.position, v_old.normal });
	}

	for (const glm::uvec3& f_old : mesh.triangles) {
		// initialize face vertices
		int f_idx = result.faces.size();
		result.faces.push_back({ f_old });

		for (int i = 0; i < 3; i++) {
			int v_idx_a = f_old[i];
			int v_idx_b = f_old[(i + 1) % 3];
			VertexEx& v_a = result.vertices[v_idx_a];
			VertexEx& v_b = result.vertices[v_idx_b];

			// initialize edge vertices
			int e_curr_idx = result.edges.size();
			for (int e_idx : v_a.edges) {
				EdgeEx& e = result.edges[e_idx];
				if (e.vertices[0] == v_idx_b || e.vertices[1] == v_idx_b) {
					e_curr_idx = e_idx;
					break;
				}
			}

			if (e_curr_idx == result.edges.size()) {
				EdgeEx new_e = { {v_idx_a, v_idx_b} };

				result.edges.push_back(new_e);

				// initialize vertex edges
				v_a.edges.push_back(e_curr_idx);
				v_b.edges.push_back(e_curr_idx);
			}

			// initialize face edges
			result.faces[f_idx].edges.push_back(e_curr_idx);

			// initialize edge faces
			EdgeEx& e_curr = result.edges[e_curr_idx];
			e_curr.faces.push_back(f_idx);
		}
	}

	return result;
}

int MeshEx::otherVertex(int e_idx, int v_idx) const {
	const EdgeEx& e = edges[e_idx];

	if (e.vertices[0] == v_idx)
		return e.vertices[1];
	else
		return e.vertices[0];
}

int MeshEx::otherEdge(int v_idx, int f_idx, int e_idx) const {
	const FaceEx& f = faces[f_idx];

	for (int e_new_idx : f.edges) {
		if (e_new_idx == e_idx)
			continue;

		const EdgeEx& e_new = edges[e_new_idx];
		if (e_new.vertices[0] == v_idx || e_new.vertices[1] == v_idx)
			return e_new_idx;
	}
}

int MeshEx::otherFace(int e_idx, int f_idx) const {
	const EdgeEx& e = edges[e_idx];

	if (e.faces[0] == f_idx)
		return e.faces[1];
	else
		return e.faces[0];
}

int MeshEx::commonEdgeOfFaces(int f_a_idx, int f_b_idx) const {
	for (int e_a_idx : faces[f_a_idx].edges) {
		for (int e_b_idx : faces[f_b_idx].edges) {
			if (e_a_idx == e_b_idx)
				return e_a_idx;
		}
	}
	return -1;
}

glm::vec3 MeshEx::vertexToVertex(int v_src_idx, int v_dst_idx) const {
	return this->vertices[v_dst_idx].position - this->vertices[v_src_idx].position;
}

double MeshEx::angleBetweenEdges(int e_a_idx, int e_b_idx) const {
	glm::uvec2 vs_a = this->edges[e_a_idx].vertices;
	glm::uvec2 vs_b = this->edges[e_b_idx].vertices;

	for (int i = 0; i < 4; i++) {
		if (vs_a[0] == vs_b[0])
			break;

		if (i % 2 == 0)
			vs_a = { vs_a[1], vs_a[0] };
		else
			vs_b = { vs_b[1], vs_b[0] };
	}

	const glm::vec3 vector_a = vertexToVertex(vs_a[0], vs_a[1]);
	const glm::vec3 vector_b = vertexToVertex(vs_b[0], vs_b[1]);

	return glm::acos(glm::dot(vector_a, vector_b) / (glm::length(vector_a) * glm::length(vector_b)));
}

// Obtained from: https://gamedev.stackexchange.com/questions/60630/how-do-i-find-the-circumcenter-of-a-triangle-in-3d
glm::vec3 MeshEx::circumcircleCenter(int f_idx) const {
	const FaceEx& f = faces[f_idx];

	glm::vec3 a = vertices[f.vertices[0]].position;
	glm::vec3 ab = vertexToVertex(f.vertices[0], f.vertices[1]);
	glm::vec3 ac = vertexToVertex(f.vertices[0], f.vertices[2]);
	glm::vec3 abXac = glm::cross(ab, ac);

	return a + (glm::cross(abXac, ab) * glm::dot(ac, ac) + glm::cross(ac, abXac) * glm::dot(ab, ab)) / (2 * glm::dot(abXac, abXac));
}

glm::vec3 MeshEx::centerOfMass(int f_idx) const {
	const FaceEx& f = faces[f_idx];

	return (vertices[f.vertices[0]].position + vertices[f.vertices[1]].position + vertices[f.vertices[2]].position) * (1.0f / 3.0f);
}

double MeshEx::defectAroundVertex(int v_idx) const {
	const VertexEx& v = vertices[v_idx];
	int v_degree = v.edges.size();

	int e_curr_idx = v.edges[0];
	int f_curr_idx = edges[e_curr_idx].faces[0];

	double curvature = 0.0;
	for (int i = 0; i < v_degree; i++) {
		int e_next_idx = otherEdge(v_idx, f_curr_idx, e_curr_idx);
		curvature += angleBetweenEdges(e_curr_idx, e_next_idx);

		e_curr_idx = e_next_idx;
		f_curr_idx = otherFace(e_curr_idx, f_curr_idx);
	}

	return curvature - 2.0 * glm::pi<double>();
}