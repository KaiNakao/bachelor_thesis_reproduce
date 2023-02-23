#include "output_distribution.hpp"

void output_distribution(const std::vector<std::vector<int>> &cny,
                         const std::vector<double> &coor,
                         const std::vector<std::vector<double>> &gauss_point,
                         const double &young,
                         const std::vector<double> &displacement_cart,
                         const std::vector<double> &displacemnet_curv) {
    const int nelement = cny.size();

    // output solution of cartesian elements
    for (int ie = 0; ie < nelement; ie++) {
        auto node_id = cny.at(ie);

        // x at nodes
        std::vector<double> xnode(3);
        for (int inode = 0; inode < 3; inode++) {
            xnode.at(inode) = coor.at(node_id.at(inode));
        }

        // ux at nodes
        std::vector<double> ux_node(3);
        for (int inode = 0; inode < 3; inode++) {
            ux_node.at(inode) = displacement_cart.at(node_id.at(inode));
        }
        std::vector<double> dvec = ux_node;

        std::string filename_u =
            "output/dist_displacement/cart/" + std::to_string(ie) + ".csv";
        std::string filename_sig =
            "output/dist_stress/cart/" + std::to_string(ie) + ".csv";
        std::ofstream ofs_u(filename_u);
        std::ofstream ofs_sig(filename_sig);

        ofs_u << std::setprecision(10);
        ofs_sig << std::setprecision(10);

        // Gaussian quadrature
        for (int ipoint = 0; ipoint < gauss_point.size(); ipoint++) {
            double r = gauss_point.at(ipoint).at(0);

            // shape function
            std::vector<double> nvec(3);
            nvec.at(0) = r * (r - 1.) / 2.;
            nvec.at(1) = r * (r + 1.) / 2.;
            nvec.at(2) = -(r + 1.) * (r - 1.);

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

            double x = 0.;
            for (int inode = 0; inode < 3; inode++) {
                x += xnode.at(inode) * nvec.at(inode);
            }

            // elastic properties
            std::vector<std::vector<double>> dmat = {{young}};

            double ux = 0.;
            for (int inode = 0; inode < 3; inode++) {
                ux += ux_node.at(inode) * nvec.at(inode);
            }

            auto sig = matvec(dmat, matvec(bmat, dvec));

            ofs_u << x << "," << ux << std::endl;
            ofs_sig << x << "," << sig.at(0) << std::endl;
        }
    }

    // read optimized parameter from text
    std::vector<std::vector<double>> alpha_vec_global(cny.size());
    std::string record;
    std::ifstream ifs("./tmp/alpha_vec_global.csv");
    for (int ie = 0; ie < cny.size(); ie++) {
        getline(ifs, record, '\n');
        std::istringstream iss(record);
        while (getline(iss, record, ',')) {
            auto v = std::stod(record);
            alpha_vec_global.at(ie).push_back(v);
        }
    }

    // output solution of curvilinear elements
    for (int ie = 0; ie < nelement; ie++) {
        // id of each node in the element
        auto node_id = cny.at(ie);

        // x at nodes
        std::vector<double> xnode(3);
        for (int inode = 0; inode < 3; inode++) {
            xnode.at(inode) = coor.at(node_id.at(inode));
        }

        const std::vector<double> alpha_vec_element = alpha_vec_global.at(ie);

        auto x_ = [&](double p) {
            std::vector<double> pvec(alpha_vec_element.size());
            for (int i = 0; i < alpha_vec_element.size(); i++) {
                pvec.at(i) = pow(p, i + 2);
            }
            return xnode.at(2) +
                   (xnode.at(1) - xnode.at(0)) *
                       (p + inner_product(pvec, alpha_vec_element));
        };

        auto dxdp_ = [&](double p) {
            std::vector<double> pvec(alpha_vec_element.size());
            for (int i = 0; i < alpha_vec_element.size(); i++) {
                pvec.at(i) = (i + 2) * pow(p, i + 1);
            }
            return (xnode.at(1) - xnode.at(0)) *
                   (1. + inner_product(pvec, alpha_vec_element));
        };

        auto d2xdp2_ = [&](double p) {
            std::vector<double> pvec(alpha_vec_element.size());
            for (int i = 0; i < alpha_vec_element.size(); i++) {
                pvec.at(i) = (i + 1) * (i + 2) * pow(p, i);
            }
            return (xnode.at(1) - xnode.at(0)) *
                   inner_product(pvec, alpha_vec_element);
        };

        auto p_ = [&](double x, double p0) {
            auto p = p0;
            double err = 1.;
            while (err > pow(10., -8)) {
                p -= (x_(p) - x) / dxdp_(p);
                err = fabs(x_(p) - x);
            }
            return p;
        };

        // ux at nodes
        std::vector<double> ux_node(3);
        for (int inode = 0; inode < 3; inode++) {
            ux_node.at(inode) = displacemnet_curv.at(node_id.at(inode));
        }

        // p at nodes
        std::vector<double> pnode(3);
        for (int inode = 0; inode < 3; inode++) {
            pnode.at(inode) = p_(xnode.at(inode), 0.);
        }

        // up2 at nodes
        std::vector<double> up_node(3);
        for (int inode = 0; inode < 3; inode++) {
            up_node.at(inode) = dxdp_(pnode.at(inode)) * ux_node.at(inode);
        }
        auto dvec = up_node;

        std::string filename_u =
            "output/dist_displacement/curv/" + std::to_string(ie) + ".csv";
        std::string filename_sig =
            "output/dist_stress/curv/" + std::to_string(ie) + ".csv";
        std::ofstream ofs_u(filename_u);
        std::ofstream ofs_sig(filename_sig);

        ofs_u << std::setprecision(10);
        ofs_sig << std::setprecision(10);

        // Gaussian quadrature
        for (int ipoint = 0; ipoint < gauss_point.size(); ipoint++) {
            double r = gauss_point.at(ipoint).at(0);

            // shape function
            std::vector<double> nvec(3);
            nvec.at(0) = r * (r - 1.) / 2.;
            nvec.at(1) = r * (r + 1.) / 2.;
            nvec.at(2) = -(r + 1.) * (r - 1.);

            std::vector<std::vector<double>> nmat(1, std::vector<double>(3));
            for (int inode = 0; inode < 3; inode++) {
                nmat.at(0).at(inode) = nvec.at(inode);
            }

            // p
            double p = 0.;
            for (int inode = 0; inode < 3; inode++) {
                p += pnode.at(inode) * nvec.at(inode);
            }

            // x
            double x = x_(p);

            // dndr (dndr[i] = dn[i]/dr)
            std::vector<double> dndr = {-1. / 2. + r, 1. / 2. + r, -2. * r};

            // dpdr (dpdr = dp/dr)
            double dpdr = 0.;
            for (int i = 0; i < 3; i++) {
                dpdr += dndr.at(i) * pnode.at(i);
            }

            // dxdp (dxdp = dx/dp)
            auto dxdp = dxdp_(p);

            // drdp (drdp = dr/dp)
            double drdp = pow(dpdr, -1.);

            // dndp (dndp[i] = dn[i]/dp)
            std::vector<double> dndp(3);
            for (int i = 0; i < 3; i++) {
                dndp.at(i) += dndr.at(i) * drdp;
            }

            // calculate bmat
            std::vector<std::vector<double>> bmat(1, std::vector<double>(3));
            for (int inode = 0; inode < 3; inode++) {
                bmat.at(0).at(inode) = dndp.at(inode);
            }

            // dpdx (dpdx = dp/dx);
            auto dpdx = pow(dxdp, -1);

            // d2xdp2 (d2xdp2 = d^2x/dp^2)
            double d2xdp2 = d2xdp2_(p);

            // elastic coefficient
            auto ctensor_cart = young;
            double ctensor = pow(dpdx, 4.) * ctensor_cart;

            // dmat (c tensor by matrix)
            std::vector<std::vector<double>> dmat(1, std::vector<double>(1));
            dmat.at(0).at(0) = ctensor;

            // calculate Christoffel symbols
            // gtensor[i][j][k] = G^k_ij
            double gtensor = d2xdp2 * dpdx;

            // gmat (g tensor by matrix)
            std::vector<std::vector<double>> gmat(1, std::vector<double>(1));
            gmat.at(0).at(0) = gtensor;

            auto gnmat = matmat(gmat, nmat);
            for (int i = 0; i < 1; i++) {
                for (int j = 0; j < 3; j++) {
                    bmat.at(i).at(j) -= gnmat.at(i).at(j);
                }
            }

            double u = 0.;
            for (int inode = 0; inode < 3; inode++) {
                u += up_node.at(inode) * nvec.at(inode);
            }
            u *= dpdx;

            auto sig = matvec(dmat, matvec(bmat, dvec));
            ofs_u << x << "," << u << std::endl;
            ofs_sig << x << "," << sig.at(0) * dxdp * dxdp << std::endl;
        }
    }
}