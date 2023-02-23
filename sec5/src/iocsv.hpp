#pragma once
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

std::vector<std::vector<int>> read_cny(const std::string &cny_path);

std::vector<std::vector<double>> read_coor(const std::string &coor_path);

std::vector<std::vector<double>> read_gauss_point(
    const std::string &gauss_point_path);

std::vector<double> read_alpha_vec(const std::string &path);

void output_displacement(const std::string &displacement_path,
                         const std::vector<std::vector<double>> &displacement,
                         const std::vector<std::vector<double>> &coor);
