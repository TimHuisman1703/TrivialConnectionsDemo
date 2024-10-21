#include "solve.h"

std::vector<std::pair<std::vector<int>, int>> getCycles(const MeshEx& mesh_ex, const std::vector<std::vector<int>>& noncon_cycles, const std::vector<int>& noncon_ks) {
    std::vector<std::pair<std::vector<int>, int>> cycles;

    // Add noncontractible cycles
    for (int i = 0; i < noncon_cycles.size(); i++)
        cycles.push_back({ noncon_cycles[i], noncon_ks[i] });

    // Add contractible cycles
    for (const VertexEx& v : mesh_ex.vertices) {
        std::vector<int> cycle;
        for (int e_idx : v.edges)
            cycle.push_back(e_idx);
        cycles.push_back({ cycle, v.k - 1 });
    }
    
    return cycles;
}

std::vector<double> calculateAdjustmentAngles(const MeshEx& mesh_ex, std::vector<std::pair<std::vector<int>, int>> cycles) {
    std::vector<Eigen::Triplet<double>> A_entries;
    std::vector<double> b_entries;
    for (int cycle_idx = 0; cycle_idx < cycles.size(); cycle_idx++) {
        const std::pair<std::vector<int>, int>& cycle_data = cycles[cycle_idx];
        const std::vector<int>& cycle_edges = cycle_data.first;
        const int k = cycle_data.second;

        for (int i = 0; i < cycle_edges.size(); i++) {
            int e_a_idx = cycle_edges[i];
            int e_b_idx = cycle_edges[(i + 1) % cycle_edges.size()];

            int f_to_idx = mesh_ex.commonFaceOfEdges(e_a_idx, e_b_idx);
            int f_from_idx = mesh_ex.otherFace(e_a_idx, f_to_idx);
            double value = f_from_idx < f_to_idx ? 1.0 : -1.0;
            A_entries.push_back(Eigen::Triplet<double>(cycle_idx, e_a_idx, value));
        }
        double defect = mesh_ex.angleOnPath(cycle_edges);

        b_entries.push_back(2 * glm::pi<double>() * k - defect);
    }

    Eigen::SparseMatrix<double> A(cycles.size(), mesh_ex.edges.size());
    A.setFromTriplets(A_entries.begin(), A_entries.end());
    Eigen::VectorXd b(cycles.size());
    for (int i = 0; i < cycles.size(); i++)
        b[i] = b_entries[i];

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A * A.transpose());
    Eigen::VectorXd u = solver.solve(b);
    Eigen::VectorXd x = A.transpose() * u;
    
    std::vector<double> adjustment_angles;
    for (int i = 0; i < x.size(); i++)
        adjustment_angles.push_back(x[i]);
    return adjustment_angles;
}