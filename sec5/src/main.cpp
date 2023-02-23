#include <cmath>

#include "calc_cart.hpp"
#include "calc_curv.hpp"
#include "iocsv.hpp"
#include "output_distribution.hpp"

int main() {
    // connectivity
    auto cny = read_cny("model/cny.csv");

    // coordinate
    auto coor = read_coor("model/coor.csv");

    // Gaussian quadrature
    auto gauss_point = read_gauss_point("gauss_point4.csv");

    // body force
    std::vector<double> fbody = {0., -15.};

    // solve in cartesian coordinates
    auto displacement_cart = calc_cart(cny, coor, gauss_point, fbody);
    output_displacement("output/displacement_node_cart.csv", displacement_cart,
                        coor);

    output_distribution_cart(cny, coor, displacement_cart,
                             "output/displacement_dist_cart.csv",
                             "output/stress_dist_cart.csv");

    // optimized parameter for coordinate transform
    auto alpha_vec = read_alpha_vec("alpha_vec_optimized_slope.csv");

    // solve in curvilinear coordinates
    auto displacement_curv =
        calc_curv(cny, coor, gauss_point, fbody, alpha_vec);
    output_displacement("output/displacement_node_curv.csv", displacement_curv,
                        coor);
    output_distribution_curv(cny, coor, displacement_curv, alpha_vec,
                             "output/displacement_dist_curv.csv",
                             "output/stress_dist_curv.csv");
}
