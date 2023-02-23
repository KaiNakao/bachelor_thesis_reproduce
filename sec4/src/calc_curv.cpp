#include "calc_curv.hpp"

std::vector<double> calc_curv(
    const std::vector<std::vector<int>> &cny,
    const std::vector<double> &coor_cart,
    const std::vector<double> &displacement,
    const std::vector<std::vector<double>> &gauss_point,
    std::function<std::vector<double>(const double &)> fbody_func,
    const double &young) {
    int nnode = coor_cart.size();

    // optimize parameter for coordinate transform
    auto alpha_vec_global =
        optimize_alpha_vec(cny, coor_cart, gauss_point, displacement);

    // initialize stiffness matrix and force vector
    std::vector<std::vector<double>> kglobal(nnode, std::vector<double>(nnode));
    std::vector<double> bglobal(nnode);

    compose_matrix_vector_curv(cny, coor_cart, fbody_func, alpha_vec_global,
                               young, gauss_point, kglobal, bglobal);

    return solve_equation(coor_cart, kglobal, bglobal);
}