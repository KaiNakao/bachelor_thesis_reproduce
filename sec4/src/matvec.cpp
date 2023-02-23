#include "matvec.hpp"

std::vector<std::vector<double>> matinv(std::vector<std::vector<double>> mat) {
    int n = mat.size();
    for (int i = 0; i < n; i++) {
        if (mat.at(i).size() != n) {
            std::cout << "size error" << std::endl;
            return {{}};
        }
    }
    std::vector<std::vector<double>> inv(n, std::vector<double>(n));
    std::vector<std::vector<double>> sweep(n, std::vector<double>(2 * n));

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sweep.at(i).at(j) = mat.at(i).at(j);
        }
        sweep.at(i).at(n + i) = 1.;
    }

    for (int k = 0; k < n; k++) {
        double max = fabs(sweep.at(k).at(k));
        int max_i = k;

        for (int i = k + 1; i < n; i++) {
            if (fabs(sweep.at(i).at(k)) > max) {
                max = fabs(sweep.at(i).at(k));
                max_i = i;
            }
        }

        if (fabs(sweep.at(max_i).at(k)) <= pow(10., -8.)) {
            std::cout << "Singular Matrix" << std::endl;
            return {{}};
        }

        if (k != max_i) {
            for (int j = 0; j < n * 2; j++) {
                double tmp = sweep.at(max_i).at(j);
                sweep.at(max_i).at(j) = sweep.at(k).at(j);
                sweep.at(k).at(j) = tmp;
            }
        }

        double a = 1. / sweep.at(k).at(k);

        for (int j = 0; j < n * 2; j++) {
            sweep.at(k).at(j) *= a;
        }

        for (int i = 0; i < n; i++) {
            if (i == k) {
                continue;
            }

            a = -sweep.at(i).at(k);

            for (int j = 0; j < n * 2; j++) {
                sweep.at(i).at(j) += sweep.at(k).at(j) * a;
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            inv.at(i).at(j) = sweep.at(i).at(n + j);
        }
    }

    return inv;
}

std::vector<double> matvec(const std::vector<std::vector<double>> &mat,
                           const std::vector<double> &vec) {
    int m = mat.size();
    int n = vec.size();
    for (int i = 0; i < m; i++) {
        if (mat.at(i).size() != n) {
            std::cout << "size error" << std::endl;
            return {};
        }
    }
    std::vector<double> ret(m);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            ret.at(i) += mat.at(i).at(j) * vec.at(j);
        }
    }

    return ret;
}

std::vector<std::vector<double>> matmat(
    const std::vector<std::vector<double>> &mat1,
    const std::vector<std::vector<double>> &mat2) {
    int m = mat1.size();
    int n = mat2.size();
    int l = mat2.at(0).size();
    for (int i = 0; i < m; i++) {
        if (mat1.at(i).size() != n) {
            std::cout << "size error" << std::endl;
            return {{}};
        }
    }
    for (int i = 0; i < n; i++) {
        if (mat2.at(i).size() != l) {
            std::cout << "size error" << std::endl;
            return {{}};
        }
    }

    std::vector<std::vector<double>> ret(m, std::vector<double>(l));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < l; j++) {
            for (int k = 0; k < n; k++) {
                ret.at(i).at(j) += mat1.at(i).at(k) * mat2.at(k).at(j);
            }
        }
    }
    return ret;
}

std::vector<std::vector<double>> transpose(
    const std::vector<std::vector<double>> &mat) {
    int m = mat.size();
    int n = mat.at(0).size();
    std::vector<std::vector<double>> ret(n, std::vector<double>(m));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            ret.at(i).at(j) = mat.at(j).at(i);
        }
    }
    return ret;
}

double matdet(const std::vector<std::vector<double>> &mat) {
    int n = mat.size();
    for (int i = 0; i < n; i++) {
        if (mat.at(i).size() != n) {
            std::cout << "size error" << std::endl;
            return 0;
        }
    }

    if (n == 2) {
        return mat.at(0).at(0) * mat.at(1).at(1) -
               mat.at(0).at(1) * mat.at(1).at(0);
    } else {
        std::cout << "size error" << std::endl;
        return 0;
    }
}

double inner_product(const std::vector<double> &vec1,
                     const std::vector<double> &vec2) {
    int n = vec1.size();
    if (vec2.size() != n) {
        std::cout << "size error" << std::endl;
        return 0;
    }
    double ret = 0.;
    for (int i = 0; i < n; i++) {
        ret += vec1.at(i) * vec2.at(i);
    }
    return ret;
}

double norm2(const std::vector<double> &vec) {
    return sqrt(inner_product(vec, vec));
}

std::vector<double> cg_solver(const std::vector<std::vector<double>> &mat,
                              const std::vector<double> &vec) {
    int n = mat.size();
    for (int i = 0; i < n; i++) {
        if (mat.at(i).size() != n) {
            std::cout << "size error" << std::endl;
            return {};
        }
    }
    if (vec.size() != n) {
        std::cout << "size error" << std::endl;
        return {};
    }

    double eps = pow(10., -10.);
    std::vector<double> x(n);
    std::vector<double> ax = matvec(mat, x);
    std::vector<double> r(n), p(n);
    for (int i = 0; i < n; i++) {
        r.at(i) = vec.at(i) - ax.at(i);
        p.at(i) = r.at(i);
    }
    double err = norm2(r) / norm2(vec);

    int cnt = 0;
    while (err > eps) {
        std::vector<double> ap = matvec(mat, p);
        double alpha = inner_product(p, r) / inner_product(p, ap);
        for (int i = 0; i < n; i++) {
            x.at(i) += alpha * p.at(i);
            r.at(i) -= alpha * ap.at(i);
        }
        double beta = inner_product(r, ap) / inner_product(p, ap);
        for (int i = 0; i < n; i++) {
            p.at(i) = r.at(i) - beta * p.at(i);
        }
        err = norm2(r) / norm2(vec);
        cnt++;
        // std::cout << err << std::endl;
    }
    // std::cout << "nloop: " << cnt << std::endl;
    return x;
}
