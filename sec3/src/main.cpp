#include <cmath>

#include "calc_cart.hpp"
#include "calc_curv.hpp"

int main() {
    const std::string input_dir("model/div5/");
    const std::string output_dir("output/div5/");

    // connectivity
    auto cny = read_cny(input_dir + "cny.csv");

    // coordinate
    auto coor = read_coor(input_dir + "coor.csv");

    // Gaussian quadrature
    auto gauss_point = read_gauss_point("gauss_point4.csv");

    // body force
    std::vector<double> fbody = {0., -15.};

    // solve in cartesion coordinates
    auto displacement_cart = calc_cart(cny, coor, gauss_point, fbody);
    output_displacement(output_dir + "displacement_node_cart.csv",
                        displacement_cart, coor);

    // solve in curvilinear coorinates
    auto displacement_curv = calc_curv(cny, coor, gauss_point, fbody);
    output_displacement(output_dir + "displacement_node_curv.csv",
                        displacement_curv, coor);
}
