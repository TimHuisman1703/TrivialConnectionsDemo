#pragma once
#include "framework/mesh.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/constants.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <optional>
#include <span>
#include <vector>
#include "mesh_ex.h";

std::vector<std::pair<std::vector<int>, int>> getCycles(const MeshEx& mesh_ex, const std::vector<std::vector<int>>& noncon_cycles, const std::vector<int>& noncon_ks);
std::vector<double> calculateAdjustmentAngles(const MeshEx& mesh_ex, std::vector<std::pair<std::vector<int>, int>> cycles);