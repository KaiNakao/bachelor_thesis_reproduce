#pragma once
#include <functional>
#include <iostream>
#include <set>
#include <vector>

#include "compose_matrix_vector_curv.hpp"
#include "optimize_alpha.hpp"
#include "iocsv.hpp"
#include "matvec.hpp"
#include "solve_equation.hpp"

std::vector<double> calc_curv(
    const std::vector<std::vector<int>> &cny, const std::vector<double> &coor_cart,
    const std::vector<double> &displacement,
    const std::vector<std::vector<double>> &gauss_point,
    std::function<std::vector<double>(const double &)> fbody_func,
    const double &young);
