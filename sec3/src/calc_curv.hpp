#pragma once
#include <functional>
#pragma once
#include <iostream>
#include <vector>

#include "compose_matrix_vector_curv.hpp"
#include "iocsv.hpp"
#include "matvec.hpp"
#include "solve_equation.hpp"

std::vector<std::vector<double>> calc_curv(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &fbody);