#include "calc_curv.hpp"

std::vector<std::vector<double>> calc_curv(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &fbody) {
    int nnode = coor.size();

    // initialize stiffness matrix and force vector
    std::vector<std::vector<double>> kglobal(2 * nnode,
                                             std::vector<double>(2 * nnode));
    std::vector<double> bglobal(2 * nnode);

    compose_matrix_vector_curv(cny, coor, gauss_point, fbody, kglobal, bglobal);
    auto ret = solve_equation(coor, kglobal, bglobal);
    return ret;
}