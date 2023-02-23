#include "compose_matrix_vector_curv.hpp"

void compose_matrix_vector_curv(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &fbody_cart,
    std::vector<std::vector<double>> &kglobal, std::vector<double> &bglobal) {
    auto p_ = [&](std::vector<double> x) {
        std::vector<double> p(2);
        p.at(0) = sqrt(x.at(0) + 1.) + sqrt(x.at(1) + 1.);
        p.at(1) = sqrt(x.at(0) + 1.) - sqrt(x.at(1) + 1.);
        return p;
    };

    // dxdp[i][j] = dxi / dpj
    auto dxdp_ = [&](std::vector<double> p) {
        return std::vector<std::vector<double>>{
            {(p.at(0) + p.at(1)) / 2., (p.at(0) + p.at(1)) / 2.},
            {(p.at(0) - p.at(1)) / 2., -(p.at(0) - p.at(1)) / 2.}};
    };

    // d2xdp2[i][j][k] = d^2xi / (dpj dpk)
    auto d2xdp2_ = [&]() {
        return std::vector<std::vector<std::vector<double>>>{
            {{1. / 2., 1. / 2.}, {1. / 2., 1. / 2.}},
            {{1. / 2., -1. / 2.}, {-1. / 2., 1. / 2.}}};
    };

    int nelement = cny.size();
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

        // physical property
        double young = 5000, nyu = 0.30;
        double lam = young * nyu / (1. + nyu) / (1. - 2. * nyu);
        double mu = young / 2. / (1. + nyu);

        // elastic coefficient cart
        std::vector<std::vector<std::vector<std::vector<double>>>> ctensor_cart(
            2, std::vector<std::vector<std::vector<double>>>(
                   2, std::vector<std::vector<double>>(
                          2, std::vector<double>(2))));
        ctensor_cart = {
            {{{2. * mu + lam, 0.}, {0., lam}}, {{0., mu}, {mu, 0.}}},
            {{{0., mu}, {mu, 0.}}, {{lam, 0.}, {0., 2. * mu + lam}}}};

        // element stiffness matrix curv
        std::vector<std::vector<double>> kelement(6, std::vector<double>(6));
        // element force vector curv
        std::vector<double> belement(6);

        // Gaussian quadrature on triangle
        for (int ipoint = 0; ipoint < gauss_point.size(); ipoint++) {
            double r1 = gauss_point.at(ipoint).at(0);
            double r2 = gauss_point.at(ipoint).at(1);
            double w = gauss_point.at(ipoint).at(2);

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
            auto d2xdp2 = d2xdp2_();

            // calculate Christoffel symbols
            // gtensor[i][j][k] = G^k_ij
            std::vector<std::vector<std::vector<double>>> gtensor(
                2, std::vector<std::vector<double>>(2, std::vector<double>(2)));
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

            // eleastic tensor curv
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

            // body force curv
            std::vector<double> fbody(2);
            for (int i = 0; i < 2; i++) {
                for (int ic = 0; ic < 2; ic++) {
                    fbody.at(i) += dpdx.at(i).at(ic) * fbody_cart.at(ic);
                }
            }

            // add to kelement, belement
            auto ktmp = matmat(transpose(bmat), matmat(dmat, bmat));
            auto btmp = matvec(transpose(nmat), fbody);
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    kelement.at(i).at(j) += ktmp.at(i).at(j) *
                                            fabs(matdet(dxdp)) *
                                            fabs(matdet(dpdr)) * w;
                }
                belement.at(i) +=
                    btmp.at(i) * fabs(matdet(dxdp)) * fabs(matdet(dpdr)) * w;
            }
        }

        // transform to cartesian coordinate
        std::vector<std::vector<double>> kelement_cart(6,
                                                       std::vector<double>(6));
        std::vector<double> belement_cart(6);
        for (int inode = 0; inode < 3; inode++) {
            auto dxdp_i = dxdp_(pnode.at(inode));
            for (int jnode = 0; jnode < 3; jnode++) {
                auto dxdp_j = dxdp_(pnode.at(jnode));
                std::vector<std::vector<double>> kblock(2,
                                                        std::vector<double>(2));
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        kblock.at(i).at(j) =
                            kelement.at(2 * inode + i).at(2 * jnode + j);
                    }
                }
                auto kblock_cart =
                    matmat(dxdp_i, matmat(kblock, transpose(dxdp_j)));
                for (int i = 0; i < 2; i++) {
                    for (int j = 0; j < 2; j++) {
                        kelement_cart.at(2 * inode + i).at(2 * jnode + j) =
                            kblock_cart.at(i).at(j);
                    }
                }
            }
            std::vector<double> bblock(2);
            for (int i = 0; i < 2; i++) {
                bblock.at(i) = belement.at(2 * inode + i);
            }
            auto bblock_cart = matvec(dxdp_i, bblock);
            for (int i = 0; i < 2; i++) {
                belement_cart.at(2 * inode + i) = bblock_cart.at(i);
            }
        }

        // assembling
        for (int i = 0; i < 2; i++) {
            for (int inode = 0; inode < 3; inode++) {
                for (int j = 0; j < 2; j++) {
                    for (int jnode = 0; jnode < 3; jnode++) {
                        kglobal.at(2 * node_id.at(inode) + i)
                            .at(2 * node_id.at(jnode) + j) +=
                            kelement_cart.at(2 * inode + i).at(2 * jnode + j);
                    }
                }
                bglobal.at(2 * node_id.at(inode) + i) +=
                    belement_cart.at(2 * inode + i);
            }
        }
    }
    return;
}