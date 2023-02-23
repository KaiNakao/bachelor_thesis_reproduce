#include "calc_cart.hpp"

std::vector<std::vector<double>> calc_cart(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &fbody) {
    int nnode = coor.size();

    // initialize stiffness matrix and force vector
    std::vector<std::vector<double>> kglobal(2 * nnode,
                                             std::vector<double>(2 * nnode));
    std::vector<double> bglobal(2 * nnode);

    // compose stiffness matrix and force vector
    compose_matrix_vector_cart(cny, coor, gauss_point, fbody, kglobal, bglobal);

    auto ret = solve_equation(coor, kglobal, bglobal);

    return ret;
}
