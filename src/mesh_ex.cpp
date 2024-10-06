#include "mesh_ex.h"

MeshEx MeshEx::fromMesh(const Mesh mesh)
{
	MeshEx result;

	// Initialize vertex positions + normals
	for (const Vertex& v_old : mesh.vertices) {
		result.vertices.push_back({ v_old.position, v_old.normal });
	}

	for (const glm::uvec3& f_old : mesh.triangles) {
		// Initialize face vertices
		int f_idx = result.faces.size();
		result.faces.push_back({ f_old });

		for (int i = 0; i < 3; i++) {
			int v_idx_a = f_old[i];
			int v_idx_b = f_old[(i + 1) % 3];
			VertexEx& v_a = result.vertices[v_idx_a];
			VertexEx& v_b = result.vertices[v_idx_b];

			// Initialize edge vertices
			int e_curr_idx = result.edges.size();
			for (int e_idx : v_a.edges) {
				EdgeEx& e = result.edges[e_idx];
				if (e.vertices[0] == v_idx_b || e.vertices[1] == v_idx_b) {
					e_curr_idx = e_idx;
					break;
				}
			}

			bool need_new_edge = e_curr_idx == result.edges.size();
			if (need_new_edge) {
				EdgeEx new_e = { {v_idx_a, v_idx_b} };

				result.edges.push_back(new_e);

				// Initialize vertex edges
				v_a.edges.push_back(e_curr_idx);
				v_b.edges.push_back(e_curr_idx);

			}

			// Initialize edge faces
			EdgeEx& e_curr = result.edges[e_curr_idx];
			e_curr.faces[1 - need_new_edge] = f_idx;
		}
	}

	// Initialize face edges
	for (int f_idx = 0; f_idx < result.faces.size(); f_idx++) {
		FaceEx& f = result.faces[f_idx];
		for (int i = 0; i < 3; i++) {
			int v_a_idx = f.vertices[i];
			int v_b_idx = f.vertices[(i + 1) % 3];
			int e_idx = result.commonEdgeOfVertices(v_a_idx, v_b_idx);
			f.edges.push_back(e_idx);
		}
	}

	// Re-order vertex edges (counter-clockwise)
	for (int v_idx = 0; v_idx < result.vertices.size(); v_idx++) {
		VertexEx& v = result.vertices[v_idx];
		int v_degree = v.edges.size();

		int e_curr_idx = v.edges[0];
		bool facing_outwards = result.edges[0].vertices[0] == v_idx;
		int f_curr_idx = result.edges[e_curr_idx].faces[facing_outwards];

		std::vector<int> edges_ordered = {};
		for (int i = 0; i < v_degree; i++) {
			edges_ordered.push_back(e_curr_idx);

			e_curr_idx = result.otherEdge(v_idx, f_curr_idx, e_curr_idx);;
			f_curr_idx = result.otherFace(e_curr_idx, f_curr_idx);
		}

		v.edges = edges_ordered;
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

int MeshEx::commonVertexOfEdges(int e_a_idx, int e_b_idx) const {
	for (int a = 0; a < 2; a++) {
		int v_a_idx = edges[e_a_idx].vertices[a];
		for (int b = 0; b < 2; b++) {
			int v_b_idx = edges[e_b_idx].vertices[b];
			if (v_a_idx == v_b_idx)
				return v_a_idx;
		}
	}
	return -1;
}

int MeshEx::commonEdgeOfVertices(int v_a_idx, int v_b_idx) const {
	for (int e_a_idx : vertices[v_a_idx].edges) {
		for (int e_b_idx : vertices[v_b_idx].edges) {
			if (e_a_idx == e_b_idx)
				return e_a_idx;
		}
	}
	return -1;
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

int MeshEx::commonFaceOfEdges(int e_a_idx, int e_b_idx) const {
	for (int a = 0; a < 2; a++) {
		int f_a_idx = edges[e_a_idx].faces[a];
		for (int b = 0; b < 2; b++) {
			int f_b_idx = edges[e_b_idx].faces[b];
			if (f_a_idx == f_b_idx)
				return f_a_idx;
		}
	}
	return -1;
}

glm::vec3 MeshEx::vertexToVertex(int v_src_idx, int v_dst_idx) const {
	return vertices[v_dst_idx].position - vertices[v_src_idx].position;
}

glm::vec3 MeshEx::edgeVector(int e_idx) const {
	const EdgeEx& e = edges[e_idx];

	return vertexToVertex(e.vertices[0], e.vertices[1]);
}

double MeshEx::angleBetweenEdges(int e_a_idx, int e_b_idx) const {
	glm::vec3 vector_a = edgeVector(e_a_idx);
	glm::vec3 vector_b = edgeVector(e_b_idx) * (edgesFormChain(e_a_idx, e_b_idx) ? -1.0f : 1.0f);

	return glm::acos(glm::dot(vector_a, vector_b) / (glm::length(vector_a) * glm::length(vector_b)));
}

double MeshEx::signedAngleBetweenEdges(int e_a_idx, int e_b_idx) const {
	return angleBetweenEdges(e_a_idx, e_b_idx) * edgeOrientation(e_a_idx, e_b_idx);
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

int MeshEx::edgesFormChain(int e_a_idx, int e_b_idx) const {
	const EdgeEx& e_a = edges[e_a_idx];
	const EdgeEx& e_b = edges[e_b_idx];
	return e_a.vertices[0] == e_b.vertices[1] || e_a.vertices[1] == e_b.vertices[0];
}

int MeshEx::edgeOrientation(int e_a_idx, int e_b_idx) const {
	int f_idx = commonFaceOfEdges(e_a_idx, e_b_idx);
	const FaceEx& f = faces[f_idx];

	int e_a_loc = f.edges[0] == e_a_idx ? 0 : f.edges[1] == e_a_idx ? 1 : 2;
	int e_b_loc = f.edges[0] == e_b_idx ? 0 : f.edges[1] == e_b_idx ? 1 : 2;

	bool ccw = (3 + e_a_loc - e_b_loc) % 3 == 1;

	return 2 * ccw - 1;
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

double MeshEx::angleOnPath(std::vector<int> path) {
	double angle = 0.0f;
	for (int i = 0; i < path.size(); i++) {
		int e_a_idx = path[i];
		int e_b_idx = path[(i + 1) % path.size()];

		angle += signedAngleBetweenEdges(e_a_idx, e_b_idx);
	}

	return angle;
}