#include "optimize_alpha.hpp"
std::vector<std::vector<double>> optimize_alpha_vec(
    const std::vector<std::vector<int>> &cny, const std::vector<double> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &displacement) {
    // directory for tmp files
    system("mkdir -p ./tmp");

    std::ofstream ofs_cny("./tmp/cny.csv");
    for (int i = 0; i < cny.size(); i++) {
        for (int j = 0; j < cny.at(i).size(); j++) {
            ofs_cny << cny.at(i).at(j);
            if (j < cny.at(i).size() - 1) {
                ofs_cny << ",";
            } else {
                ofs_cny << std::endl;
            }
        }
    }

    // coordinate of nodes
    std::ofstream ofs_coor("./tmp/coor.csv");
    ofs_coor << std::setprecision(10);
    for (int i = 0; i < coor.size(); i++) {
        ofs_coor << coor.at(i) << std::endl;
    }

    // gausian quadrature
    std::ofstream ofs_gausspoint("./tmp/gauss_point.csv");
    ofs_gausspoint << std::setprecision(10);
    for (int i = 0; i < gauss_point.size(); i++) {
        ofs_gausspoint << gauss_point.at(i).at(0) << ","
                       << gauss_point.at(i).at(1) << std::endl;
    }

    // solution in cartesian coordinates
    std::ofstream ofs_displacement("./tmp/displacement.csv");
    ofs_displacement << std::setprecision(10);
    for (int i = 0; i < displacement.size(); i++) {
        ofs_displacement << displacement.at(i) << std::endl;
    }

    // execute python script for parameter optimization
    system("python3 src/optimize_alpha_vec.py");

    // read optimized parameter from text
    std::vector<std::vector<double>> alpha_vec_global(cny.size());
    std::string record;
    std::ifstream ifs("./tmp/alpha_vec_global.csv");
    for (int ie = 0; ie < cny.size(); ie++) {
        getline(ifs, record, '\n');
        std::istringstream iss(record);
        while (getline(iss, record, ',')) {
            auto v = std::stod(record);
            alpha_vec_global.at(ie).push_back(v);
        }
    }

    // system("rm -r ./tmp");
    return alpha_vec_global;
}
