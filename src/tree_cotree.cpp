#pragma once
#include "tree_cotree.h"

std::vector<int> treeCotreeDecompose(const MeshEx& mesh_ex) {
	std::vector<int> tree_assignment(mesh_ex.edges.size());

	// Create cotree
	std::queue<int> f_queue;
	f_queue.push(0);
	std::vector<bool> f_seen(mesh_ex.faces.size());
	f_seen[0] = true;

	while (!f_queue.empty()) {
		int f_curr_idx = f_queue.front();
		f_queue.pop();

		for (int e_idx : mesh_ex.faces[f_curr_idx].edges) {
			int f_next_idx = mesh_ex.otherFace(e_idx, f_curr_idx);

			if (tree_assignment[e_idx] == 0 && !f_seen[f_next_idx]) {
				f_seen[f_next_idx] = true;
				f_queue.push(f_next_idx);
				tree_assignment[e_idx] = -1;
			}
		}
	}

	// Create tree
	std::queue<int> v_queue;
	v_queue.push(0);
	std::vector<bool> v_seen(mesh_ex.vertices.size());
	v_seen[0] = true;

	while (!v_queue.empty()) {
		int v_curr_idx = v_queue.front();
		v_queue.pop();

		for (int e_idx : mesh_ex.vertices[v_curr_idx].edges) {
			int v_next_idx = mesh_ex.otherVertex(e_idx, v_curr_idx);

			if (tree_assignment[e_idx] == 0 && !v_seen[v_next_idx]) {
				v_seen[v_next_idx] = true;
				v_queue.push(v_next_idx);
				tree_assignment[e_idx] = 1;
			}
		}
	}

	return tree_assignment;
}

std::vector<std::vector<int>> findNoncontractibleCycles(const MeshEx& mesh_ex, std::vector<int> tree_assignment) {
	std::vector<std::vector<int>> paths;

	for (int e_idx = 0; e_idx < mesh_ex.edges.size(); e_idx++) {
		if (tree_assignment[e_idx] != 0)
			continue;

		const EdgeEx& e = mesh_ex.edges[e_idx];
		int f_start_idx = e.faces[0];
		int f_end_idx = e.faces[1];

		std::queue<int> f_queue;
		f_queue.push(f_start_idx);
		std::vector<int> came_from;
		for (int f_idx = 0; f_idx < mesh_ex.faces.size(); f_idx++)
			came_from.push_back(-1);
		came_from[f_start_idx] = f_end_idx;

		while (!f_queue.empty()) {
			int f_curr_idx = f_queue.front();
			f_queue.pop();

			if (f_curr_idx == f_end_idx)
				break;

			for (int e_idx : mesh_ex.faces[f_curr_idx].edges) {
				int f_next_idx = mesh_ex.otherFace(e_idx, f_curr_idx);

				if (tree_assignment[e_idx] == -1 && came_from[f_next_idx] == -1 && came_from[f_curr_idx] != f_next_idx) {
					came_from[f_next_idx] = f_curr_idx;
					f_queue.push(f_next_idx);
				}
			}
		}

		std::vector<int> path;
		int f_prev_idx = f_start_idx;
		int f_curr_idx = f_end_idx;
		while (true) {
			int e_idx = mesh_ex.commonEdgeOfFaces(f_prev_idx, f_curr_idx);
			path.push_back(e_idx);

			f_prev_idx = f_curr_idx;
			f_curr_idx = came_from[f_curr_idx];
			if (f_curr_idx == f_end_idx)
				break;
		}

		paths.push_back(path);
	}

	return paths;
}