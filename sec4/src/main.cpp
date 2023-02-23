#include "calc_cart.hpp"
#include "calc_curv.hpp"
#include "iocsv.hpp"
#include "output_distribution.hpp"

int main() {
    const std::string cny_path = "model/cny.csv";
    const std::string coor_path = "model/coor.csv";
    const std::string gauss_point_path = "gauss_point.csv";

    // connectivity
    auto cny = read_cny(cny_path);

    // coordinate
    auto coor = read_coor(coor_path);

    // Gaussian quadrature
    auto gauss_point = read_gauss_point(gauss_point_path);

    // distribution of body force
    auto fbody_ = [&](const double &x) {
        return std::vector<double>{(x - 25.) * (x - 75.) / 50.};
    };

    // elastic coefficient
    double young = 5000.;

    // solve in cartesian coordinates
    auto displacement_cart = calc_cart(cny, coor, gauss_point, fbody_, young);
    output_displacement("output/displacement_cart.csv", displacement_cart,
                        coor);

    // solve in curvilinear coordinates
    auto displacement_curv =
        calc_curv(cny, coor, displacement_cart, gauss_point, fbody_, young);
    output_displacement("output/displacement_curv.csv", displacement_curv,
                        coor);

    output_distribution(cny, coor, gauss_point, young, displacement_cart,
                        displacement_curv);
}
