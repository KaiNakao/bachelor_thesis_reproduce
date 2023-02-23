#pragma once
#include <iomanip>

#include "matvec.hpp"

void output_distribution(const std::vector<std::vector<int>> &cny,
                         const std::vector<double> &coor,
                         const std::vector<std::vector<double>> &gauss_point,
                         const double &young,
                         const std::vector<double> &displacement_cart,
                         const std::vector<double> &displacemnet_curv);