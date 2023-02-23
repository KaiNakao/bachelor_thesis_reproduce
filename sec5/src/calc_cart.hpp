#pragma once
#include <functional>
#include <iostream>
#include <vector>

#include "compose_matrix_vector_cart.hpp"
#include "iocsv.hpp"
#include "matvec.hpp"
#include "solve_equation.hpp"

std::vector<std::vector<double>> calc_cart(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &fbody);