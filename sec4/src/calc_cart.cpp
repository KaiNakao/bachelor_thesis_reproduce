#include "calc_cart.hpp"

std::vector<double> calc_cart(
    const std::vector<std::vector<int>> &cny, const std::vector<double> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    std::function<std::vector<double>(const double &)> fbody_func,
    const double &young) {
    int nnode = coor.size();

    // initialize stiffness matrix and force vector
    std::vector<std::vector<double>> kglobal(nnode, std::vector<double>(nnode));
    std::vector<double> bglobal(nnode);

    compose_matrix_vector_cart(cny, coor, fbody_func, young, gauss_point,
                               kglobal, bglobal);

    return solve_equation(coor, kglobal, bglobal);
}
