import numpy as np
from scipy.optimize import minimize


def lglobal(alpha_vec):
    xmin = 0.
    xmax = 8.
    xc = (xmin + xmax) / 2.
    dx = xmax - xmin
    ymin = 0.
    ymax = 5.
    yc = (ymin + ymax) / 2.
    dy = ymax - ymin

    alphax = alpha_vec[:9].copy()
    alphay = alpha_vec[9:].copy()
    alphax[0] += 1
    alphay[1] += 1

    def x_(p):
        x = np.zeros(2)
        pvec = np.array([p[0], p[1], p[0]**2, p[0]*p[1], p[1] **
                        2, p[0]**3, p[0]**2 * p[1], p[0] * p[1]**2, p[1]**3])
        xr = np.dot(alphax, pvec)
        yr = np.dot(alphay, pvec)
        x[0] = xc + xr * dx
        x[1] = yc + yr * dy
        return x

    def dxdp_(p):
        dxdp = np.zeros((2, 2))
        dpvec = np.array([[1., 0., 2. * p[0], p[1], 0., 3. * p[0]**2, 2. * p[0] * p[1], p[1]**2, 0.],
                          [0., 1., 0., p[0], 2. * p[1], 0., p[0]**2, 2. * p[0] * p[1], 3. * p[1]**2]])
        for i in range(2):
            dxdp[0][i] = np.dot(alphax, dpvec[i]) * dx
            dxdp[1][i] = np.dot(alphay, dpvec[i]) * dy
        return dxdp

    def d2xdp2_(p):
        d2xdp2 = np.zeros((2, 2, 2))
        ddpvec = np.array([[[0., 0., 2., 0., 0., 6. * p[0], 2. * p[1], 0., 0.],
                            [0., 0., 0., 1., 0., 0., 2. * p[0], 2. * p[1], 0.]],
                           [[0., 0., 0., 1., 0., 0., 2. * p[0], 2. * p[1], 0.],
                            [0., 0., 0., 0., 2., 0., 0., 2. * p[0], 6. * p[1]]]])
        for i in range(2):
            for j in range(2):
                d2xdp2[0][i][j] = np.dot(alphax, ddpvec[i][j]) * dx
                d2xdp2[1][i][j] = np.dot(alphay, ddpvec[i][j]) * dy
        return d2xdp2

    def p_(x):
        p0 = np.array([(x[0] - xc) / dx, (x[1] - yc) / dy])
        p = p0
        err = 1.
        cnt = 0
        while (err > pow(10, -8) and cnt < 10):
            delta = x_(p) - x
            dp = np.dot(np.linalg.inv(dxdp_(p)), delta)
            p -= dp
            err = np.linalg.norm(dp)
            cnt += 1
        return [p, cnt]

    # def x_(p):
    #     x = np.zeros(2)
    #     x[0] = pow((p[0] + p[1]) / 2., 2) - 1.
    #     x[1] = pow((p[0] - p[1]) / 2., 2) - 1.
    #     return x

    # def p_(x):
    #     p = np.zeros(2)
    #     p[0] = np.sqrt(x[0] + 1.) + np.sqrt(x[1] + 1.)
    #     p[1] = np.sqrt(x[0] + 1.) - np.sqrt(x[1] + 1.)
    #     return p

    # def dxdp_(p): # dxdp[i][j] = dxi / dpj
    #     return [[(p[0] + p[1])/2., (p[0] + p[1])/2.], [(p[0] - p[1])/2., -(p[0] - p[1])/2.]]

    # def d2xdp2_(p): # d2xdp2[i][j][k] = d^2xi / (dpj dpk)
    #     return [[[1./2., 1./2.], [1./2., 1./2.]], [[1./2., -1./2.], [-1./2., 1./2.]]]

    nelement = len(cny)
    lglobal = 0.
    for ie in range(nelement):
        # std::cout << "element #" << ie << std::endl
        # id of each node in the element
        node_id = cny[ie]

        # (x, y) at nodes
        xnode = np.zeros((3, 2))
        for inode in range(3):
            xnode[inode] = coor[node_id[inode]]

        # (p, q) at nodes
        pnode = np.zeros((3, 2))
        for inode in range(3):
            pnode[inode], cnt = p_(xnode[inode])
            if cnt == 10:
                return 100
        pvec = np.zeros(6)
        for inode in range(3):
            pvec[2 * inode + 0] = pnode[inode][0]
            pvec[2 * inode + 1] = pnode[inode][1]

        # ux at nodes
        uxnode = np.zeros((3, 2))
        for inode in range(3):
            uxnode[inode] = displacement[node_id[inode]]
        # up at nodes
        upnodes = np.zeros((3, 2))
        for inode in range(3):
            dxdp = dxdp_(pnode[inode])
            ux = uxnode[inode]
            up = np.zeros(2)
            for i in range(2):
                for j in range(2):
                    up[i] += dxdp[j][i] * ux[j]
            upnodes[inode] = up
        dvec = np.zeros(6)
        for inode in range(3):
            dvec[2 * inode + 0] = upnodes[inode][0]
            dvec[2 * inode + 1] = upnodes[inode][1]

        # physical property
        young = young_global[ie]
        nyu = 0.30
        lam = young * nyu / (1. + nyu) / (1. - 2. * nyu)
        mu = young / 2. / (1. + nyu)

        # elastic coefficient
        ctensor_cart = np.array([[[[2. * mu + lam, 0.], [0., lam]],
                                  [[0., mu], [mu, 0.]]],
                                 [[[0., mu], [mu, 0.]],
                                  [[lam, 0.], [0., 2. * mu + lam]]]])

        lelement = 0.
        # Gaussian quadrature on triangle
        for ipoint in range(len(gauss_point)):
            r1 = gauss_point[ipoint][0]
            r2 = gauss_point[ipoint][1]
            w = gauss_point[ipoint][2]

            # shape function
            nvec = np.array([1. - r1 - r2, r1, r2])

            nmat = np.zeros((2, 6))
            for inode in range(3):
                nmat[0][2 * inode] = nvec[inode]
                nmat[1][2 * inode + 1] = nvec[inode]

            # dndr (dndr[i][j] = dn[i]/dr[j])
            dndr = [[-1., -1.], [1., 0.], [0., 1.]]

            # dpdr (dpdr[i][j] = dp[i]/dr[j])
            dpdr = np.zeros((2, 2))
            for i in range(2):
                for j in range(2):
                    for n in range(3):
                        dpdr[i][j] += pnode[n][i] * dndr[n][j]

            # drdp (drdp[i][j] = dr[i]/dp[j])
            drdp = np.linalg.inv(dpdr)

            # dndp (dndp[i][j] = dn[i]/dp[j])
            dndp = np.zeros((3, 2))
            for i in range(3):
                for j in range(2):
                    for k in range(2):
                        dndp[i][j] += dndr[i][k] * drdp[k][j]

            # calculate bmat
            bmat = np.zeros((4, 6))
            for inode in range(3):
                bmat[0][inode * 2 + 0] = dndp[inode][0]
                bmat[1][inode * 2 + 0] = dndp[inode][1] / 2.
                bmat[1][inode * 2 + 1] = dndp[inode][0] / 2.
                bmat[2][inode * 2 + 0] = dndp[inode][1] / 2.
                bmat[2][inode * 2 + 1] = dndp[inode][0] / 2.
                bmat[3][inode * 2 + 1] = dndp[inode][1]

            # derivative of coordinate transformation
            p = np.zeros(2)
            for inode in range(3):
                p[0] += pnode[inode][0] * nvec[inode]
                p[1] += pnode[inode][1] * nvec[inode]
            dxdp = dxdp_(p)
            dpdx = np.linalg.inv(dxdp)
            d2xdp2 = d2xdp2_(p)

            # calculate Christoffel symbols
            gtensor = np.zeros((2, 2, 2))  # gtensor[i][j][k] = G^k_ij
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for m in range(2):
                            gtensor[i][j][k] += d2xdp2[m][i][j] * dpdx[k][m]
            gmat = np.zeros((4, 2))
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        gmat[2 * i + j][k] = gtensor[i][j][k]

            gnmat = np.dot(gmat, nmat)
            for i in range(4):
                for j in range(6):
                    bmat[i][j] -= gnmat[i][j]

            ctensor = np.zeros((2, 2, 2, 2))
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for l in range(2):
                            for ic in range(2):
                                for jc in range(2):
                                    for kc in range(2):
                                        for lc in range(2):
                                            ctensor[i][j][k][l] += dpdx[i][ic] * dpdx[j][jc] * \
                                                dpdx[k][kc] * dpdx[l][lc] * \
                                                ctensor_cart[ic][jc][kc][lc]

            dmat = np.zeros((4, 4))
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for l in range(2):
                            dmat[2 * i + j][2 * k + l] = ctensor[i][j][k][l]

            # body force
            fbody = np.zeros(2)
            for i in range(2):
                for ic in range(2):
                    fbody[i] += dpdx[i][ic] * fbody_cart[ic]

            # up
            up = np.dot(nmat, dvec)

            # strain
            epsvec = np.dot(bmat, dvec)
            epstensor = np.zeros((2, 2))
            for i in range(2):
                for j in range(2):
                    epstensor[i][j] = epsvec[2 * i + j]

            l1 = 0.
            l2 = 0.
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        for l in range(2):
                            l1 += epstensor[i][j] * \
                                ctensor[i][j][k][l] * epstensor[k][l] / 2.
                l2 += up[i] * fbody[i]
            lelement += (l1 - l2) * abs(np.linalg.det(dxdp)) * \
                abs(np.linalg.det(dpdr)) * w
        lglobal += lelement
    print(lglobal, np.linalg.norm(alpha_vec)**2)
    lam = 5.
    lglobal += lam * np.linalg.norm(alpha_vec)**2
    return lglobal


if __name__ == "__main__":
    cny = np.loadtxt("model/slope_model/cny.csv",
                     skiprows=1, delimiter=",", dtype=int)[:, 1:]
    coor = np.loadtxt("model/slope_model/coor.csv",
                      skiprows=1, delimiter=",")[:, 1:]
    young_global = np.loadtxt(
        "model/slope_model/young.csv", skiprows=1, delimiter=",")[:, 1]
    displacement = np.loadtxt(
        "output/slope/displacement_node_cart.csv", skiprows=1, delimiter=",")[:, 2:]
    gauss_point = np.loadtxt("gauss_point4.csv", skiprows=1, delimiter=",")
    fbody_cart = np.array([0., -15.])
    alpha_init = np.zeros(18)
    res = minimize(lglobal, alpha_init, method="nelder-mead",
                   options={"maxiter": 10000})
    print("result: ", lglobal(alpha_init), " -> ", lglobal(res.x))
    print(res)
    np.savetxt("alpha_vec_optimized_slope.csv", res.x, delimiter=",")
