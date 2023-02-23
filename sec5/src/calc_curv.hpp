#pragma once
#include <functional>
#include <iostream>
#include <vector>

#include "compose_matrix_vector_curv.hpp"
#include "iocsv.hpp"
#include "matvec.hpp"
#include "output_distribution.hpp"
#include "solve_equation.hpp"

std::vector<std::vector<double>> calc_curv(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &fbody, const std::vector<double> &alpha_vec);