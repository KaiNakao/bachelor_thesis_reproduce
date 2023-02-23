#pragma once
#pragma once
#include <functional>
#include <iomanip>

#include "iocsv.hpp"
#include "matvec.hpp"

std::vector<std::vector<double>> optimize_alpha_vec(
    const std::vector<std::vector<int>> &cny, const std::vector<double> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &displacement);