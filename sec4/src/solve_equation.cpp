#include "solve_equation.hpp"

std::vector<double> solve_equation(const std::vector<double> &coor,
                                   std::vector<std::vector<double>> &kglobal,
                                   std::vector<double> &bglobal) {
    const int nnode = coor.size();
    // Reorder the equations
    int cnt0 = 0, cnt1 = 0;
    std::vector<int> trans(nnode);
    for (int inode = 0; inode < nnode; inode++) {
        double x = coor.at(inode);
        // boundary
        if (x == 0) {
            trans.at(inode) = nnode - 1 - cnt0;
            cnt0++;
        } else {
            trans.at(inode) = cnt1;
            cnt1++;
        }
    }

    std::vector<std::vector<double>> k_reordered(nnode,
                                                 std::vector<double>(nnode));
    std::vector<double> b_reordered(nnode);
    for (int i = 0; i < nnode; i++) {
        for (int j = 0; j < nnode; j++) {
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
    std::vector<std::vector<double>> kab(m, std::vector<double>(nnode - m));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < nnode - m; j++) {
            kab.at(i).at(j) = k_reordered.at(i).at(m + j);
        }
    }
    std::vector<double> fa(m);
    for (int i = 0; i < m; i++) {
        fa.at(i) = b_reordered.at(i);
    }

    // Given value at boundary
    std::vector<double> ub(nnode - m, 0);

    std::vector<double> kaaub = matvec(kab, ub);
    std::vector<double> vec(m);
    for (int i = 0; i < m; i++) {
        vec.at(i) = fa.at(i) - kaaub.at(i);
    }

    // Solve the equation (unknown value vector ua)
    std::vector<double> ua = cg_solver(kaa, vec);

    // Organize the solution
    std::vector<double> tmp(nnode);
    for (int i = 0; i < m; i++) {
        tmp.at(i) = ua.at(i);
    }
    for (int i = 0; i < nnode - m; i++) {
        tmp.at(m + i) = ub.at(i);
    }
    std::vector<double> result(nnode);
    for (int inode = 0; inode < nnode; inode++) {
        result.at(inode) = tmp.at(trans.at(inode));
    }

    return result;
}