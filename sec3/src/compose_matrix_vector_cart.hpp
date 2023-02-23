#pragma once
#include <fstream>
#pragma once
#include <functional>
#include <sstream>
#include <vector>

#include "matvec.hpp"

void compose_matrix_vector_cart(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &fbody, std::vector<std::vector<double>> &kglobal,
    std::vector<double> &bglobal);
