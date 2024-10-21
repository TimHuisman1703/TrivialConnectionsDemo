#pragma once
#include "framework/mesh.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/constants.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <optional>
#include <span>
#include <vector>


struct VertexEx {
	glm::vec3 position;
	glm::vec3 normal;
	int k = 0;

	std::vector<int> edges = {};

	[[nodiscard]] constexpr bool operator==(const VertexEx&) const noexcept = default;
};

struct EdgeEx {
	glm::uvec2 vertices;

	glm::uvec2 faces = { -1, -1 };

	[[nodiscard]] constexpr bool operator==(const EdgeEx&) const noexcept = default;
};

struct FaceEx {
	glm::uvec3 vertices;
	glm::vec3 normal;

	std::vector<int> edges = {};

	[[nodiscard]] constexpr bool operator==(const FaceEx&) const noexcept = default;
};

struct MeshEx {
	std::vector<VertexEx> vertices = {};
	std::vector<EdgeEx> edges = {};
	std::vector<FaceEx> faces = {};

	static MeshEx fromMesh(const Mesh);

	int otherVertex(int e_idx, int v_idx) const;
	int otherEdge(int v_idx, int f_idx, int e_idx) const;
	int otherFace(int e_idx, int f_idx) const;

	int commonVertexOfEdges(int e_a_idx, int e_b_idx) const;
	int commonEdgeOfVertices(int v_a_idx, int v_b_idx) const;
	int commonEdgeOfFaces(int f_a_idx, int f_b_idx) const;
	int commonFaceOfEdges(int e_a_idx, int e_b_idx) const;

	glm::vec3 vertexToVertex(int v_src_idx, int v_dst_idx) const;
	glm::vec3 edgeVector(int e_idx) const;

	glm::vec3 circumcircleCenter(int f_idx) const;
	glm::vec3 centerOfMass(int f_idx) const;

	int edgePointToEdge(int e_a_idx, int e_b_idx) const;
	int edgeOrientationInTriangle(int e_a_idx, int e_b_idx) const;

	double turnAngleBetweenEdges(int e_a_idx, int e_b_idx) const;
	double innerAngleBetweenEdges(int e_a_idx, int e_b_idx) const;
	double angleOnPath(std::vector<int> path) const;
	double angleOnPathAdjusted(std::vector<int> path, std::vector<double> adjustment_angles) const;
};