#include "iocsv.hpp"

#include <stdio.h>

std::vector<std::vector<int>> read_cny(const std::string &cny_path) {
    std::vector<std::vector<int>> ret;
    std::string record;
    std::ifstream ifs(cny_path);
    getline(ifs, record, '\n');
    while (getline(ifs, record, '\n')) {
        std::vector<int> row(3);
        std::istringstream iss(record);
        getline(iss, record, ',');
        for (int i = 0; i < row.size(); i++) {
            getline(iss, record, ',');
            row.at(i) = std::stoi(record);
        }
        ret.push_back(row);
    }
    return ret;
}

std::vector<std::vector<double>> read_coor(const std::string &coor_path) {
    std::vector<std::vector<double>> ret;
    std::string record;
    std::ifstream ifs(coor_path);
    getline(ifs, record, '\n');
    while (getline(ifs, record, '\n')) {
        std::vector<double> row(2);
        std::istringstream iss(record);
        getline(iss, record, ',');
        for (int i = 0; i < row.size(); i++) {
            getline(iss, record, ',');
            row.at(i) = std::stod(record);
        }
        ret.push_back(row);
    }
    return ret;
}

std::vector<std::vector<double>> read_gauss_point(
    const std::string &gauss_point_path) {
    std::vector<std::vector<double>> ret;
    std::string record;
    std::ifstream ifs(gauss_point_path);
    getline(ifs, record, '\n');
    while (getline(ifs, record, '\n')) {
        std::vector<double> row(3);
        std::istringstream iss(record);
        for (int i = 0; i < row.size(); i++) {
            getline(iss, record, ',');
            row.at(i) = std::stod(record);
        }
        ret.push_back(row);
    }
    return ret;
}

void output_displacement(const std::string &displacement_path,
                         const std::vector<std::vector<double>> &displacement,
                         const std::vector<std::vector<double>> &coor) {
    std::ofstream ofs(displacement_path);
    ofs << "x,y,ux,uy" << std::endl;
    ofs << std::setprecision(10);
    for (int i = 0; i < displacement.size(); i++) {
        double x = coor.at(i).at(0);
        double y = coor.at(i).at(1);
        double ux = displacement.at(i).at(0);
        double uy = displacement.at(i).at(1);
        ofs << x << "," << y << "," << ux << "," << uy << std::endl;
    }
    return;
}
