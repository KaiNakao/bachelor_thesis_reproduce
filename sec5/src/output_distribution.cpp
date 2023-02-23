#include "output_distribution.hpp"

void output_distribution_cart(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &displacement,
    const std::string &displacement_path, const std::string &stress_path) {
    std::vector<std::vector<double>> plt_point = {
        {0., 0.},    {0.25, 0.},   {0.5, 0.},   {0.75, 0.},   {1., 0.},
        {0., 0.25},  {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0., 0.5},
        {0.25, 0.5}, {0.5, 0.5},   {0., 0.75},  {0.25, 0.75}, {0., 1.}};
    std::ofstream ofs_displacement(displacement_path);
    std::ofstream ofs_stress(stress_path);
    ofs_displacement << std::setprecision(10);
    ofs_stress << std::setprecision(10);
    int nelement = cny.size();
    // loop over elements
    for (int ie = 0; ie < nelement; ie++) {
        std::cout << "output element #" << ie << std::endl;
        // id of each node in the element
        auto node_id = cny.at(ie);
        // (x, y) at nodes
        std::vector<std::vector<double>> xnode(3, std::vector<double>(2));
        for (int inode = 0; inode < 3; inode++) {
            xnode.at(inode) = coor.at(node_id.at(inode));
        }
        std::vector<double> xvec(6);
        for (int inode = 0; inode < 3; inode++) {
            xvec.at(2 * inode + 0) = xnode.at(inode).at(0);
            xvec.at(2 * inode + 1) = xnode.at(inode).at(1);
        }

        // u at nodes
        std::vector<std::vector<double>> unode(3, std::vector<double>(2));
        for (int inode = 0; inode < 3; inode++) {
            unode.at(inode) = displacement.at(node_id.at(inode));
        }
        std::vector<double> dvec(6);
        for (int inode = 0; inode < 3; inode++) {
            dvec.at(2 * inode + 0) = unode.at(inode).at(0);
            dvec.at(2 * inode + 1) = unode.at(inode).at(1);
        }
        // physical property
        double young = 5000., nyu = 0.30;
        double lam = young * nyu / (1. + nyu) / (1. - 2. * nyu);
        double mu = young / 2. / (1. + nyu);

        // elastic coefficient
        std::vector<std::vector<std::vector<std::vector<double>>>> ctensor(
            2, std::vector<std::vector<std::vector<double>>>(
                   2, std::vector<std::vector<double>>(
                          2, std::vector<double>(2))));
        ctensor = {{{{2. * mu + lam, 0.}, {0., lam}}, {{0., mu}, {mu, 0.}}},
                   {{{0., mu}, {mu, 0.}}, {{lam, 0.}, {0., 2. * mu + lam}}}};

        std::vector<std::vector<double>> dmat(4, std::vector<double>(4));
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    for (int l = 0; l < 2; l++) {
                        dmat.at(2 * i + j).at(2 * k + l) =
                            ctensor.at(i).at(j).at(k).at(l);
                    }
                }
            }
        }

        // Gaussian quadrature on triangle
        for (int ipoint = 0; ipoint < plt_point.size(); ipoint++) {
            double r1 = plt_point.at(ipoint).at(0);
            double r2 = plt_point.at(ipoint).at(1);

            // shape function
            std::vector<double> nvec = {1. - r1 - r2, r1, r2};
            std::vector<std::vector<double>> nmat(2, std::vector<double>(6));
            for (int inode = 0; inode < 3; inode++) {
                nmat.at(0).at(2 * inode) = nvec.at(inode);
                nmat.at(1).at(2 * inode + 1) = nvec.at(inode);
            }

            // dndr (dndr[i][j] = dn[i]/dr[j])
            std::vector<std::vector<double>> dndr = {
                {-1., -1.}, {1., 0.}, {0., 1.}};

            // dxdr (dxdr[i][j] = dx[i]/dr[j])
            std::vector<std::vector<double>> dxdr(2, std::vector<double>(2));
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int n = 0; n < 3; n++) {
                        dxdr.at(i).at(j) +=
                            xnode.at(n).at(i) * dndr.at(n).at(j);
                    }
                }
            }

            // drdx (drdx[i][j] = dr[i]/dx[j])
            std::vector<std::vector<double>> drdx = matinv(dxdr);

            // dndx (dndx[i][j] = dn[i]/dx[j])
            std::vector<std::vector<double>> dndx(3, std::vector<double>(2));
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        dndx.at(i).at(j) += dndr.at(i).at(k) * drdx.at(k).at(j);
                    }
                }
            }

            // calculate bmat
            std::vector<std::vector<double>> bmat(4, std::vector<double>(6));
            for (int inode = 0; inode < 3; inode++) {
                bmat.at(0).at(inode * 2 + 0) = dndx.at(inode).at(0);
                bmat.at(1).at(inode * 2 + 0) = dndx.at(inode).at(1) / 2.;
                bmat.at(1).at(inode * 2 + 1) = dndx.at(inode).at(0) / 2.;
                bmat.at(2).at(inode * 2 + 0) = dndx.at(inode).at(1) / 2.;
                bmat.at(2).at(inode * 2 + 1) = dndx.at(inode).at(0) / 2.;
                bmat.at(3).at(inode * 2 + 1) = dndx.at(inode).at(1);
            }

            auto x = matvec(nmat, xvec);
            auto u = matvec(nmat, dvec);
            auto sigvec = matvec(dmat, matvec(bmat, dvec));

            ofs_displacement << x.at(0) << " " << x.at(1) << " " << u.at(0)
                             << " " << u.at(1) << std::endl;
            ofs_stress << x.at(0) << " " << x.at(1) << " " << sigvec.at(0)
                       << " " << sigvec.at(1) << " " << sigvec.at(2) << " "
                       << sigvec.at(3) << std::endl;
        }
    }
    return;
}

void output_distribution_curv(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &displacement,
    const std::vector<double> &alpha_vec, const std::string &displacement_path,
    const std::string &stress_path) {
    std::vector<std::vector<double>> plt_point = {
        {0., 0.},    {0.25, 0.},   {0.5, 0.},   {0.75, 0.},   {1., 0.},
        {0., 0.25},  {0.25, 0.25}, {0.5, 0.25}, {0.75, 0.25}, {0., 0.5},
        {0.25, 0.5}, {0.5, 0.5},   {0., 0.75},  {0.25, 0.75}, {0., 1.}};
    std::ofstream ofs_displacement(displacement_path);
    std::ofstream ofs_stress(stress_path);
    ofs_displacement << std::setprecision(10);
    ofs_stress << std::setprecision(10);
    int nelement = cny.size();

    double xmin = 0.;
    double xmax = 8.;
    double xc = (xmin + xmax) / 2.;
    double dx = xmax - xmin;
    double ymin = 0.;
    double ymax = 5.;
    double yc = (ymin + ymax) / 2.;
    double dy = ymax - ymin;

    std::vector<double> alphax(9);
    std::vector<double> alphay(9);
    for (int i = 0; i < 9; i++) {
        alphax.at(i) = alpha_vec.at(i);
        alphay.at(i) = alpha_vec.at(9 + i);
    }
    alphax.at(0) += 1.;
    alphay.at(1) += 1.;

    auto x_ = [&](std::vector<double> p) {
        std::vector<double> x(2);
        std::vector<double> pvec = {p.at(0),
                                    p.at(1),
                                    pow(p.at(0), 2),
                                    p.at(0) * p.at(1),
                                    pow(p.at(1), 2),
                                    pow(p.at(0), 3),
                                    pow(p.at(0), 2) * p.at(1),
                                    p.at(0) * pow(p.at(1), 2),
                                    pow(p.at(1), 3)};
        double xr = inner_product(alphax, pvec);
        double yr = inner_product(alphay, pvec);
        x[0] = xc + xr * dx;
        x[1] = yc + yr * dy;
        return x;
    };

    // dxdp[i][j] = dxi / dpj
    auto dxdp_ = [&](std::vector<double> p) {
        std::vector<std::vector<double>> dxdp(2, std::vector<double>(2));
        std::vector<std::vector<double>> dpvec = {
            {1., 0., 2. * p.at(0), p.at(1), 0., 3. * pow(p.at(0), 2),
             2. * p.at(0) * p.at(1), pow(p.at(1), 2), 0.},
            {0., 1., 0., p.at(0), 2. * p.at(1), 0., pow(p.at(0), 2),
             2. * p.at(0) * p.at(1), 3. * pow(p.at(1), 2)}};
        for (int i = 0; i < 2; i++) {
            dxdp.at(0).at(i) = inner_product(alphax, dpvec.at(i)) * dx;
            dxdp.at(1).at(i) = inner_product(alphay, dpvec.at(i)) * dy;
        }
        return dxdp;
    };

    // d2xdp2[i][j][k] = d^2xi / (dpj dpk)
    auto d2xdp2_ = [&](std::vector<double> p) {
        std::vector<std::vector<std::vector<double>>> d2xdp2(
            2, std::vector<std::vector<double>>(2, std::vector<double>(2)));
        std::vector<std::vector<std::vector<double>>> ddpvec = {
            {{0., 0., 2., 0., 0., 6. * p.at(0), 2. * p.at(1), 0., 0.},
             {0., 0., 0., 1., 0., 0., 2. * p.at(0), 2. * p.at(1), 0.}},
            {{0., 0., 0., 1., 0., 0., 2. * p.at(0), 2. * p.at(1), 0.},
             {0., 0., 0., 0., 2., 0., 0., 2. * p.at(0), 6. * p.at(1)}}};
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                d2xdp2.at(0).at(i).at(j) =
                    inner_product(alphax, ddpvec.at(i).at(j)) * dx;
                d2xdp2.at(1).at(i).at(j) =
                    inner_product(alphay, ddpvec.at(i).at(j)) * dy;
            }
        }
        return d2xdp2;
    };

    // solve x(p) = x for p
    auto p_ = [&](std::vector<double> x) {
        std::vector<double> p0 = {(x.at(0) - xc) / dx, (x.at(1) - yc) / dy};
        auto p = p0;
        double err = 1.;
        int cnt = 0;
        while (err > pow(10, -8) || cnt < 10) {
            auto xp = x_(p);
            std::vector<double> delta(2);
            for (int i = 0; i < 2; i++) {
                delta.at(i) = xp.at(i) - x.at(i);
            }
            auto dp = matvec(matinv(dxdp_(p)), delta);
            for (int i = 0; i < 2; i++) {
                p.at(i) -= dp.at(i);
            }
            err = norm2(dp);
            cnt++;
        }
        return p;
    };

    for (int ie = 0; ie < nelement; ie++) {
        // std::cout << "element #" << ie << std::endl;
        // id of each node in the element
        auto node_id = cny.at(ie);

        // (x, y) at nodes
        std::vector<std::vector<double>> xnode(3, std::vector<double>(2));
        for (int inode = 0; inode < 3; inode++) {
            xnode.at(inode) = coor.at(node_id.at(inode));
        }
        // (p, q) at nodes
        std::vector<std::vector<double>> pnode(3, std::vector<double>(2));
        for (int inode = 0; inode < 3; inode++) {
            pnode.at(inode) = p_(xnode.at(inode));
        }
        std::vector<double> pvec(6);
        for (int inode = 0; inode < 3; inode++) {
            pvec.at(2 * inode + 0) = pnode.at(inode).at(0);
            pvec.at(2 * inode + 1) = pnode.at(inode).at(1);
        }

        // ux at nodes
        std::vector<std::vector<double>> uxnode(3, std::vector<double>(2));
        for (int inode = 0; inode < 3; inode++) {
            uxnode.at(inode) = displacement.at(node_id.at(inode));
        }
        // up at nodes
        std::vector<std::vector<double>> upnodes(3, std::vector<double>(2));
        for (int inode = 0; inode < 3; inode++) {
            auto dxdp = dxdp_(pnode.at(inode));
            auto ux = uxnode.at(inode);
            std::vector<double> up(2);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    up.at(i) += dxdp.at(j).at(i) * ux.at(j);
                }
            }
            upnodes.at(inode) = up;
        }
        std::vector<double> dvec(6);
        for (int inode = 0; inode < 3; inode++) {
            dvec.at(2 * inode + 0) = upnodes.at(inode).at(0);
            dvec.at(2 * inode + 1) = upnodes.at(inode).at(1);
        }

        // physical property
        double young = 5000., nyu = 0.30;
        double lam = young * nyu / (1. + nyu) / (1. - 2. * nyu);
        double mu = young / 2. / (1. + nyu);

        // elastic coefficient
        std::vector<std::vector<std::vector<std::vector<double>>>> ctensor_cart(
            2, std::vector<std::vector<std::vector<double>>>(
                   2, std::vector<std::vector<double>>(
                          2, std::vector<double>(2))));
        ctensor_cart = {
            {{{2. * mu + lam, 0.}, {0., lam}}, {{0., mu}, {mu, 0.}}},
            {{{0., mu}, {mu, 0.}}, {{lam, 0.}, {0., 2. * mu + lam}}}};

        // Gaussian quadrature on triangle
        for (int ipoint = 0; ipoint < plt_point.size(); ipoint++) {
            double r1 = plt_point.at(ipoint).at(0);
            double r2 = plt_point.at(ipoint).at(1);

            // shape function
            std::vector<double> nvec = {1. - r1 - r2, r1, r2};

            std::vector<std::vector<double>> nmat(2, std::vector<double>(6));
            for (int inode = 0; inode < 3; inode++) {
                nmat.at(0).at(2 * inode) = nvec.at(inode);
                nmat.at(1).at(2 * inode + 1) = nvec.at(inode);
            }

            // dndr (dndr[i][j] = dn[i]/dr[j])
            std::vector<std::vector<double>> dndr = {
                {-1., -1.}, {1., 0.}, {0., 1.}};

            // dpdr (dpdr[i][j] = dp[i]/dr[j])
            std::vector<std::vector<double>> dpdr(2, std::vector<double>(2));
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int n = 0; n < 3; n++) {
                        dpdr.at(i).at(j) +=
                            pnode.at(n).at(i) * dndr.at(n).at(j);
                    }
                }
            }

            // drdp (drdp[i][j] = dr[i]/dp[j])
            std::vector<std::vector<double>> drdp = matinv(dpdr);

            // dndp (dndp[i][j] = dn[i]/dp[j])
            std::vector<std::vector<double>> dndp(3, std::vector<double>(2));
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        dndp.at(i).at(j) += dndr.at(i).at(k) * drdp.at(k).at(j);
                    }
                }
            }

            // calculate bmat
            std::vector<std::vector<double>> bmat(4, std::vector<double>(6));
            for (int inode = 0; inode < 3; inode++) {
                bmat.at(0).at(inode * 2 + 0) = dndp.at(inode).at(0);
                bmat.at(1).at(inode * 2 + 0) = dndp.at(inode).at(1) / 2.;
                bmat.at(1).at(inode * 2 + 1) = dndp.at(inode).at(0) / 2.;
                bmat.at(2).at(inode * 2 + 0) = dndp.at(inode).at(1) / 2.;
                bmat.at(2).at(inode * 2 + 1) = dndp.at(inode).at(0) / 2.;
                bmat.at(3).at(inode * 2 + 1) = dndp.at(inode).at(1);
            }

            // derivative of coordinate transformation
            std::vector<double> p(2);
            for (int inode = 0; inode < 3; inode++) {
                p.at(0) += pnode.at(inode).at(0) * nvec.at(inode);
                p.at(1) += pnode.at(inode).at(1) * nvec.at(inode);
            }
            auto dxdp = dxdp_(p);
            auto dpdx = matinv(dxdp);
            auto d2xdp2 = d2xdp2_(p);

            // calculate Christoffel symbols
            std::vector<std::vector<std::vector<double>>> gtensor(
                2,
                std::vector<std::vector<double>>(
                    2, std::vector<double>(2)));  // gtensor[i][j][k] = G^k_ij
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        for (int m = 0; m < 2; m++) {
                            gtensor.at(i).at(j).at(k) +=
                                d2xdp2.at(m).at(i).at(j) * dpdx.at(k).at(m);
                        }
                    }
                }
            }
            std::vector<std::vector<double>> gmat(4, std::vector<double>(2));
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        gmat.at(2 * i + j).at(k) = gtensor.at(i).at(j).at(k);
                    }
                }
            }

            auto gnmat = matmat(gmat, nmat);
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 6; j++) {
                    bmat.at(i).at(j) -= gnmat.at(i).at(j);
                }
            }

            std::vector<std::vector<std::vector<std::vector<double>>>> ctensor(
                2, std::vector<std::vector<std::vector<double>>>(
                       2, std::vector<std::vector<double>>(
                              2, std::vector<double>(2))));
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        for (int l = 0; l < 2; l++) {
                            for (int ic = 0; ic < 2; ic++) {
                                for (int jc = 0; jc < 2; jc++) {
                                    for (int kc = 0; kc < 2; kc++) {
                                        for (int lc = 0; lc < 2; lc++) {
                                            ctensor.at(i).at(j).at(k).at(l) +=
                                                dpdx.at(i).at(ic) *
                                                dpdx.at(j).at(jc) *
                                                dpdx.at(k).at(kc) *
                                                dpdx.at(l).at(lc) *
                                                ctensor_cart.at(ic)
                                                    .at(jc)
                                                    .at(kc)
                                                    .at(lc);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            std::vector<std::vector<double>> dmat(4, std::vector<double>(4));
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        for (int l = 0; l < 2; l++) {
                            dmat.at(2 * i + j).at(2 * k + l) =
                                ctensor.at(i).at(j).at(k).at(l);
                        }
                    }
                }
            }

            // up
            auto up = matvec(nmat, dvec);

            // stress
            auto sigvec = matvec(dmat, matvec(bmat, dvec));
            std::vector<std::vector<double>> sigp(2, std::vector<double>(2));
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    sigp.at(i).at(j) = sigvec.at(2 * i + j);
                }
            }

            // x
            auto x = x_(p);

            // ux
            std::vector<double> ux(2);
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    ux.at(i) += dpdx.at(j).at(i) * up.at(j);
                }
            }

            // sigx
            std::vector<std::vector<double>> sigx(2, std::vector<double>(2));
            for (int i = 0; i < 2; i++) {
                for (int j = 0; j < 2; j++) {
                    for (int k = 0; k < 2; k++) {
                        for (int l = 0; l < 2; l++) {
                            sigx.at(i).at(j) += dxdp.at(i).at(k) *
                                                dxdp.at(j).at(l) *
                                                sigp.at(k).at(l);
                        }
                    }
                }
            }

            ofs_displacement << x.at(0) << " " << x.at(1) << " " << ux.at(0)
                             << " " << ux.at(1) << std::endl;
            ofs_stress << x.at(0) << " " << x.at(1) << " " << sigx.at(0).at(0)
                       << " " << sigx.at(0).at(1) << " " << sigx.at(1).at(0)
                       << " " << sigx.at(1).at(1) << std::endl;
        }
    }
    return;
}