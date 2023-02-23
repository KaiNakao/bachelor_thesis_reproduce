#include "compose_matrix_vector_cart.hpp"

void compose_matrix_vector_cart(
    const std::vector<std::vector<int>> &cny,
    const std::vector<std::vector<double>> &coor,
    const std::vector<std::vector<double>> &gauss_point,
    const std::vector<double> &fbody, std::vector<std::vector<double>> &kglobal,
    std::vector<double> &bglobal) {
    int nelement = cny.size();

    // loop over elements
    for (int ie = 0; ie < nelement; ie++) {
        // std::cout << "element #" << ie << std::endl;
        // id of each node in the element
        auto node_id = cny.at(ie);

        // (x, y) at nodes
        std::vector<std::vector<double>> xnode(3, std::vector<double>(2));
        for (int inode = 0; inode < 3; inode++) {
            xnode.at(inode) = coor.at(node_id.at(inode));
        }

        // element stiffness matrix
        std::vector<std::vector<double>> kelement(6, std::vector<double>(6));
        // element force vector
        std::vector<double> belement(6);

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
        for (int ipoint = 0; ipoint < gauss_point.size(); ipoint++) {
            double r1 = gauss_point.at(ipoint).at(0);
            double r2 = gauss_point.at(ipoint).at(1);
            double w = gauss_point.at(ipoint).at(2);

            // shape function
            std::vector<double> nvec = {1. - r1 - r2, r1, r2};

            std::vector<std::vector<double>> nmat(2, std::vector<double>(6));
            for (int i = 0; i < 3; i++) {
                nmat.at(0).at(2 * i) = nvec.at(i);
                nmat.at(1).at(2 * i + 1) = nvec.at(i);
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

            // add to kelement, belement
            std::vector<std::vector<double>> ktmp =
                matmat(matmat(transpose(bmat), dmat), bmat);
            std::vector<double> btmp = matvec(transpose(nmat), fbody);
            for (int i = 0; i < 6; i++) {
                for (int j = 0; j < 6; j++) {
                    ktmp.at(i).at(j) *= fabs(matdet(dxdr)) / 2. * w;
                    kelement.at(i).at(j) += ktmp.at(i).at(j);
                }
                btmp.at(i) *= fabs(matdet(dxdr)) / 2. * w;
                belement.at(i) += btmp.at(i);
            }
        }

        // assembling
        for (int i = 0; i < 2; i++) {
            for (int inode = 0; inode < 3; inode++) {
                for (int j = 0; j < 2; j++) {
                    for (int jnode = 0; jnode < 3; jnode++) {
                        kglobal.at(2 * node_id.at(inode) + i)
                            .at(2 * node_id.at(jnode) + j) +=
                            kelement.at(2 * inode + i).at(2 * jnode + j);
                    }
                }
                bglobal.at(2 * node_id.at(inode) + i) +=
                    belement.at(2 * inode + i);
            }
        }
    }
}
