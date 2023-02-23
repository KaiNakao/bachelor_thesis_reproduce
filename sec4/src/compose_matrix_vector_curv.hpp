#pragma once
#include <functional>

#include "matvec.hpp"

void compose_matrix_vector_curv(
    const std::vector<std::vector<int>> &cny,
    const std::vector<double> &coor_cart,
    std::function<std::vector<double>(const double &)> fbody_func,
    const std::vector<std::vector<double>> &alpha_vec_global,
    const double &young, const std::vector<std::vector<double>> &gauss_point,
    std::vector<std::vector<double>> &kglobal, std::vector<double> &bglobal);
