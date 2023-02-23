#pragma once
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "iocsv.hpp"
#include "matvec.hpp"

void output_distribution_cart(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &displacement,
    const std::string &displacement_path, const std::string &stress_path);

void output_distribution_curv(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &displacement,
    const std::vector<double> &alpha_vec, const std::string &displacement_path,
    const std::string &stress_path);