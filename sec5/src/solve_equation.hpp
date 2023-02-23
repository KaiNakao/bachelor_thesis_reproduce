#pragma once
#include "matvec.hpp"

std::vector<std::vector<double>> solve_equation(
    const std::vector<std::vector<double>> &coor,
    std::vector<std::vector<double>> &kglobal, std::vector<double> &bglobal);