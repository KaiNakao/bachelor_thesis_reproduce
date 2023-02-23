#include "compose_matrix_vector_cart.hpp"

void compose_matrix_vector_cart(
    const std::vector<std::vector<int>> &cny, const std::vector<double> &coor,
    std::function<std::vector<double>(const double &)> fbody_func,
    const double &young, const std::vector<std::vector<double>> &gauss_point,
    std::vector<std::vector<double>> &kglobal, std::vector<double> &bglobal) {
    int nelement = cny.size();

    // loop over elements (2nd order element)
    for (int ie = 0; ie < nelement; ie++) {
        // std::cout << "element #" << ie << std::endl;
        // id of each node in the element
        auto node_id = cny.at(ie);

        // x at nodes
        std::vector<double> xnode(3);
        for (int inode = 0; inode < 3; inode++) {
            xnode.at(inode) = coor.at(node_id.at(inode));
        }

        // element stiffness matrix
        std::vector<std::vector<double>> kelement(3, std::vector<double>(3));
        // element force vector
        std::vector<double> belement(3);

        // Gaussian quadrature on triangle
        for (int ipoint = 0; ipoint < gauss_point.size(); ipoint++) {
            double r = gauss_point.at(ipoint).at(0);
            double w = gauss_point.at(ipoint).at(1);

            // shape function
            std::vector<double> nvec(3);
            nvec.at(0) = r * (r - 1.) / 2.;
            nvec.at(1) = r * (r + 1.) / 2.;
            nvec.at(2) = -(r + 1.) * (r - 1.);

            std::vector<std::vector<double>> nmat(1, std::vector<double>(3));
            for (int inode = 0; inode < 3; inode++) {
                nmat.at(0).at(inode) = nvec.at(inode);
            }

            // dndr (dndr[i] = dn[i]/dr)
            std::vector<double> dndr = {-1. / 2. + r, 1. / 2. + r, -2. * r};

            // dxdr (dxdr = dx/dr)
            double dxdr = 0;
            for (int i = 0; i < 3; i++) {
                dxdr += dndr.at(i) * xnode.at(i);
            }

            // drdx (drdx = dr/dx)
            double drdx = pow(dxdr, -1.);

            // dndx (dndx[i] = dn[i]/dx)
            std::vector<double> dndx(3);
            for (int i = 0; i < 3; i++) {
                dndx.at(i) += dndr.at(i) * drdx;
            }

            // calculate bmat
            std::vector<std::vector<double>> bmat(1, std::vector<double>(3));
            for (int inode = 0; inode < 3; inode++) {
                bmat.at(0).at(inode) = dndx.at(inode);
            }

            // body force
            double x = 0.;
            for (int inode = 0; inode < 3; inode++) {
                x += xnode.at(inode) * nvec.at(inode);
            }
            auto fbody = fbody_func(x);

            // elastic properties
            std::vector<std::vector<double>> dmat = {{young}};

            // add to kelement, belement
            std::vector<std::vector<double>> ktmp =
                matmat(matmat(transpose(bmat), dmat), bmat);
            std::vector<double> btmp = matvec(transpose(nmat), fbody);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    ktmp.at(i).at(j) *= fabs(dxdr) * w;
                    kelement.at(i).at(j) += ktmp.at(i).at(j);
                }
                btmp.at(i) *= fabs(dxdr) * w;
                belement.at(i) += btmp.at(i);
            }
        }

        // assembling
        for (int i = 0; i < 1; i++) {
            for (int inode = 0; inode < 3; inode++) {
                for (int j = 0; j < 1; j++) {
                    for (int jnode = 0; jnode < 3; jnode++) {
                        kglobal.at(node_id.at(inode) + i)
                            .at(node_id.at(jnode) + j) +=
                            kelement.at(inode + i).at(jnode + j);
                    }
                }
                bglobal.at(node_id.at(inode) + i) += belement.at(inode + i);
            }
        }
    }

    return;
}
