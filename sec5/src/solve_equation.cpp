#include "solve_equation.hpp"

std::vector<std::vector<double>> solve_equation(
    const std::vector<std::vector<double>> &coor,
    std::vector<std::vector<double>> &kglobal, std::vector<double> &bglobal) {
    int nnode = coor.size();
    // Reorder the equations
    int cnt0 = 0, cnt1 = 0;
    std::vector<int> trans(2 * nnode);
    for (int inode = 0; inode < nnode; inode++) {
        double y = coor.at(inode).at(1);
        if (fabs(y) <= pow(10., -4)) {
            trans.at(2 * inode + 0) = 2 * nnode - 1 - cnt0;
            cnt0++;
            trans.at(2 * inode + 1) = 2 * nnode - 1 - cnt0;
            cnt0++;
        } else {
            trans.at(2 * inode + 0) = cnt1;
            cnt1++;
            trans.at(2 * inode + 1) = cnt1;
            cnt1++;
        }
    }
    std::vector<std::vector<double>> k_reordered(
        2 * nnode, std::vector<double>(2 * nnode));
    std::vector<double> b_reordered(2 * nnode);
    for (int i = 0; i < 2 * nnode; i++) {
        for (int j = 0; j < 2 * nnode; j++) {
            k_reordered.at(trans.at(i)).at(trans.at(j)) = kglobal.at(i).at(j);
        }
        b_reordered.at(trans.at(i)) = bglobal.at(i);
    }

    // Define submatrix and subvector
    int m = cnt1;  // degree of freedom
    std::vector<std::vector<double>> kaa(m, std::vector<double>(m));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            kaa.at(i).at(j) = k_reordered.at(i).at(j);
        }
    }
    std::vector<std::vector<double>> kab(m, std::vector<double>(2 * nnode - m));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < 2 * nnode - m; j++) {
            kab.at(i).at(j) = k_reordered.at(i).at(m + j);
        }
    }
    std::vector<double> fa(m);
    for (int i = 0; i < m; i++) {
        fa.at(i) = b_reordered.at(i);
    }

    // Given values at boundary
    std::vector<double> ub(2 * nnode - m);

    std::vector<double> kaaub = matvec(kab, ub);
    std::vector<double> vec(m);
    for (int i = 0; i < m; i++) {
        vec.at(i) = fa.at(i) - kaaub.at(i);
    }
    // Solve the equation (unknown value vector ua)
    std::vector<double> ua = cg_solver(kaa, vec);

    // Organize the solution
    std::vector<double> tmp(2 * nnode);
    for (int i = 0; i < m; i++) {
        tmp.at(i) = ua.at(i);
    }
    for (int i = 0; i < 2 * nnode - m; i++) {
        tmp.at(m + i) = ub.at(i);
    }
    std::vector<std::vector<double>> result(nnode, std::vector<double>(2));
    for (int inode = 0; inode < nnode; inode++) {
        result.at(inode).at(0) = tmp.at(trans.at(2 * inode));
        result.at(inode).at(1) = tmp.at(trans.at(2 * inode + 1));
    }

    return result;
}