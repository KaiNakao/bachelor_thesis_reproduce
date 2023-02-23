import numpy as np
from scipy.optimize import minimize


def fbody_(x):
    # distribution of body force
    return np.array([(x - 25.) * (x - 75.)/50.])


def young_(x):
    # elastic coefficient
    return 5000.


def lelement(alpha_vec_element, xnode, ux_node, gauss_point, fbody_func, young_func):
    # calculate energy integration over one element

    def x_(p):
        pvec = np.zeros(len(alpha_vec_element))
        for i in range(len(alpha_vec_element)):
            pvec[i] = p**(i + 2)
        return xnode[2] + (xnode[1] - xnode[0]) * (p + np.dot(pvec, alpha_vec_element))

    def dxdp_(p):
        pvec = np.zeros(len(alpha_vec_element))
        for i in range(len(alpha_vec_element)):
            pvec[i] = (i + 2) * p**(i + 1)
        return (xnode[1] - xnode[0]) * (1. + np.dot(pvec, alpha_vec_element))

    def d2xdp2_(p):
        pvec = np.zeros(len(alpha_vec_element))
        for i in range(len(alpha_vec_element)):
            pvec[i] = (i + 1) * (i + 2) * p**i
        return (xnode[1] - xnode[0]) * np.dot(pvec, alpha_vec_element)

    # solve x(p) = x for p
    def p_(x, p0):
        p = p0
        err = 1.
        cnt = 0
        while (err > 10.**-8 and cnt < 10):
            p -= (x_(p) - x) / dxdp_(p)
            err = abs(x_(p) - x)
            cnt += 1
        return p

    # p at nodes
    pnode = np.array([-0.5, 0.5, 0.])
    for inode in range(3):
        pnode[inode] = p_(xnode[inode], pnode[inode])

    # up at nodes
    dvec = np.zeros(3)
    for inode in range(3):
        dxdp = dxdp_(pnode[inode])
        dvec[inode] = dxdp * ux_node[inode]

    # energy in element
    lelement = 0.

    # Gaussian quadrature on triangle
    for ipoint in range(len(gauss_point)):
        r = gauss_point[ipoint][0]
        w = gauss_point[ipoint][1]

        # shape function
        nvec = np.zeros(3)
        nvec[0] = r * (r - 1.) / 2.
        nvec[1] = r * (r + 1.) / 2.
        nvec[2] = -(r + 1.) * (r - 1.)

        nmat = np.zeros((1, 3))
        for inode in range(3):
            nmat[0][inode] = nvec[inode]

        # dndr (dndr[i] = dn[i]/dr)
        dndr = np.array([-1. / 2. + r, 1. / 2. + r, -2. * r])

        #  dpdr (dpdr = dp/dr)
        dpdr = 0.
        for inode in range(3):
            dpdr += dndr[inode] * pnode[inode]

        # drdp (drdp = dr/dp)
        drdp = dpdr**(-1.)

        # dndp (dndp[i] = dn[i]/dp)
        dndp = np.zeros(3)
        for inode in range(3):
            dndp[inode] += dndr[inode] * drdp

        # p
        p = 0.
        for inode in range(3):
            p += pnode[inode] * nvec[inode]

        # x
        x = x_(p)

        # up
        up = 0.
        for inode in range(3):
            up += nvec[inode] * dvec[inode]

        # dxdp (dxdp = dx/dp)
        dxdp = dxdp_(p)

        # dpdx (dpdx = dp/dx);
        dpdx = dxdp**(-1)

        dxdr = dxdp * dpdr

        # d2xdp2 (d2xdp2 = d^2x/dp^2)
        d2xdp2 = d2xdp2_(p)

        # calculate bmat
        bmat = np.zeros((1, 3))
        for inode in range(3):
            bmat[0][inode] = dndp[inode]

        # calculate Christoffel symbols
        # gtensor[i][j][k] = G^k_ij
        gtensor = d2xdp2 * dpdx

        # gmat (g tensor by matrix)
        gmat = np.zeros((1, 1))
        gmat[0][0] = gtensor

        gnmat = np.dot(gmat, nmat)
        bmat -= gnmat

        # elastic coefficient
        ctensor_cart = young_func(x)
        ctensor = dpdx**4. * ctensor_cart

        # dmat
        dmat = np.array([[ctensor]])

        # body force
        fbody_cart = fbody_func(x)
        fbody = np.array([dpdx * fbody_cart[0]])

        # strain
        epsvec = np.dot(bmat, dvec)

        # stress
        sigvec = np.dot(dmat, epsvec)

        # add to lelement
        lelement += (epsvec[0] * sigvec[0] / 2. -
                     fbody[0] * up) * abs(dxdr) * w
    return lelement


if __name__ == "__main__":
    coor = np.loadtxt("./tmp/coor.csv", delimiter=",")
    cny = np.loadtxt("./tmp/cny.csv", delimiter=",", dtype=int)
    gauss_point = np.loadtxt("./tmp/gauss_point.csv", delimiter=",")
    displacement = np.loadtxt("./tmp/displacement.csv", delimiter=",")

    # parameter for all elements
    alpha_vec_global = []

    # minimize energy for each element
    for ie in range(len(cny)):

        node_id = cny[ie]

        # x at nodes
        xnode = np.zeros(3)
        for inode in range(3):
            xnode[inode] = coor[node_id[inode]]

        # ux at nodes
        ux_node = np.zeros(3)
        for inode in range(3):
            ux_node[inode] = displacement[node_id[inode]]

        # optimize alpha
        def lelement_wrap(alpha_vec_element):
            le = lelement(alpha_vec_element, xnode, ux_node,
                          gauss_point, fbody_, young_)
            lam = 1.
            # for i in range(len(alpha_vec_element)):
            #     le += lam * (alpha_vec_element[i])**2
            return le

        alpha_init = np.zeros(3)
        res = minimize(lelement_wrap, alpha_init, method="BFGS")

        print(ie, res)

        alpha_vec_global.append(res.x)

    alpha_vec_global = np.array(alpha_vec_global)
    # output result
    np.savetxt("./tmp/alpha_vec_global.csv", alpha_vec_global, delimiter=",")
