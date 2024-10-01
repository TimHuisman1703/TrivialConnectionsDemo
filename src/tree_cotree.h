#pragma once
#include "mesh_ex.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/constants.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <optional>
#include <queue>
#include <span>
#include <stack>
#include <vector>


std::vector<int> treeCotreeDecompose(const MeshEx& mesh_ex);
std::vector<std::vector<int>> findNoncontractibleCycles(const MeshEx& mesh_ex, std::vector<int> tree_assignment);