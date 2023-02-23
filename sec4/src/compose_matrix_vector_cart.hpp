#pragma once
#include <functional>
#include <iomanip>

#include "matvec.hpp"

void compose_matrix_vector_cart(
    const std::vector<std::vector<int>> &cny, const std::vector<double> &coor,
    std::function<std::vector<double>(const double &)> fbody_func,
    const double &young, const std::vector<std::vector<double>> &gauss_point,
    std::vector<std::vector<double>> &kglobal, std::vector<double> &bglobal);