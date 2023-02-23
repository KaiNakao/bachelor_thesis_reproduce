#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

std::vector<std::vector<double>> matinv(std::vector<std::vector<double>> mat);

std::vector<double> matvec(const std::vector<std::vector<double>> &mat,
                           const std::vector<double> &vec);

std::vector<std::vector<double>> matmat(
    const std::vector<std::vector<double>> &mat1,
    const std::vector<std::vector<double>> &mat2);

std::vector<std::vector<double>> transpose(
    const std::vector<std::vector<double>> &mat);

double matdet(const std::vector<std::vector<double>> &mat);

double inner_product(const std::vector<double> &vec1,
                     const std::vector<double> &vec2);

double norm2(const std::vector<double> &vec);

std::vector<double> cg_solver(const std::vector<std::vector<double>> &mat,
                              const std::vector<double> &vec);
