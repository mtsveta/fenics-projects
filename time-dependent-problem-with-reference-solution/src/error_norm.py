from dolfin.cpp.mesh import cells
import numpy

__author__ = 'svetlana'

from dolfin import *
import majorant

def grad_error_norm(A, grad_u_k, grad_u_k12, grad_u_k1, grad_v_k, grad_v_k1, tau, mesh):

    int_grad_v_sq_var = tau / 3 * (inner(A * grad_v_k, grad_v_k)
                                   + inner(A * grad_v_k1, grad_v_k)
                                   + inner(A * grad_v_k1, grad_v_k1))
    int_grad_vu_var = tau / 6 * (inner(A * grad_u_k, grad_v_k)
                                 + inner(A * grad_u_k1, grad_v_k1)
                                 + 2 * inner(A * (grad_v_k1 + grad_v_k), grad_u_k12))

    int_grad_u_sq_var = tau / 15 * (2 * (inner(A * grad_u_k, grad_u_k) +
                                         inner(A * grad_u_k1, grad_u_k1))
                                    + 8 * inner(A * grad_u_k12, grad_u_k12)
                                    + 2 * inner(A * (grad_u_k1 + grad_u_k), grad_u_k12)
                                    - inner(A * grad_u_k1, grad_u_k))

    # (nabla(u - v))^2
    e_d_var = int_grad_v_sq_var - 2 * int_grad_vu_var + int_grad_u_sq_var
    e_d = assemble(e_d_var * dx(domain=mesh))

    return e_d, e_d_var

def error_norm_nd_t(k, e_incr_k, ed_k, et_k,
                    u_k1, grad_u_k, grad_u_k12, grad_u_k1, quadrature,
                    v_k1, grad_v_k, grad_v_k1,
                    A,
                    V_exact, VV_exact,
                    mesh, dim, tau, delta):

    v_k1_Ve = project(v_k1, V_exact)
    grad_v_k_Ve = project(grad_v_k, VV_exact)
    grad_v_k1_Ve = project(grad_v_k1, VV_exact)

    grad_u_k_Ve = project(grad_u_k, VV_exact)
    grad_u_k12_Ve = project(grad_u_k12, VV_exact)
    grad_u_k1_Ve = project(grad_u_k1, VV_exact)

    #e_d, e_d_var = grad_error_norm(A, grad_u_k_Ve, grad_u_k12_Ve, grad_u_k1_Ve, quadrature, V_exact, VV_exact, grad_v_k_Ve, grad_v_k1_Ve, tau, mesh)
    e_d, e_d_var = grad_error_norm(A, grad_u_k_Ve, grad_u_k12_Ve, grad_u_k1_Ve, grad_v_k, grad_v_k1, tau, mesh)
    e_t = assemble((inner(u_k1 - v_k1_Ve, u_k1 - v_k1_Ve)) * dx(domain=mesh))

    ed_k[k] = e_d
    et_k[k] = e_t

    # Calculate the error increment
    e_incr_k[k] = (2.0 - delta) * e_d + e_t

    return e_incr_k, ed_k, et_k, e_d_var

def v_norm_nd_t(k, v_norm_incr_k, vd_k, vt_k, v_k1, grad_v_k, grad_v_k1, A, V_exact, VV_exact, mesh, dim, tau, delta):

    v_k1_Ve = interpolate(v_k1, V_exact)

    grad_v_k_Ve = project(grad_v_k, VV_exact)
    grad_v_k1_Ve = project(grad_v_k1, VV_exact)

    vd_k[k] = assemble(tau / Constant(3.0) * (inner(grad_v_k_Ve, grad_v_k1_Ve) + inner(grad_v_k1_Ve, grad_v_k_Ve) + inner(grad_v_k1_Ve, grad_v_k1_Ve)) * dx(domain=mesh))
    vt_k[k] = assemble((inner(v_k1_Ve, v_k1_Ve)) * dx(domain=mesh))

    # Calculate the error increment
    v_norm_incr_k[k] = (2.0 - delta) * vd_k[k] + vt_k[k]

    return v_norm_incr_k, vd_k, vt_k



class SpaceIntegration1d():

    def __init__(self, n, dim):
        self.dim = dim
        # order
        self.n = n
        # gauss points and coefficients
        if n == 2:
            self.xi = [-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)]
            self.omega = [1.0, 1.0]
        elif n == 3:
            self.xi = [-sqrt(0.6), 0, sqrt(0.6)]
            self.omega = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
        elif n == 4:
            self.xi = [-0.861136311594953, -0.339981043584856, 0.339981043584856, 0.861136311594953]
            self.omega = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]

    def integrate_1d(self, func, diam, x1, x2):
        f_int = 0
        # gauss
        for i in xrange(self.n):
            f_int += 0.5 * diam * self.omega[i] * \
                     func(0.5 * (x1 * (1 - self.xi[i]) + x2 * (1 + self.xi[i])))
        return f_int

class SpaceIntegration2d():
    def __init__(self, n):
        # order
        self.n = n
        # gauss
        if n == 1:
            self.xi = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]
            self.omega = [1.0]
        elif n == 3:
            self.xi = [[0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [0.5, 0.0, 0.5]]
            self.omega = [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0]
        elif n == 4:
            self.xi = [[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
                       [0.6, 0.2, 0.2],
                       [0.2, 0.6, 0.2],
                       [0.2, 0.2, 0.6]]
            self.omega = [-27.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0, 25.0 / 48.0]
        elif n == 6:
            alpha_1 = 0.0597158717
            beta_1 = 0.4701420641
            alpha_2 = 0.7974269853
            beta_2 = 0.1012865073

            self.xi = [[1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0],
                       [alpha_1, beta_1, beta_1],
                       [beta_1, alpha_1, beta_1],
                       [beta_1, beta_1, alpha_1],
                       [alpha_2, beta_2, beta_2],
                       [beta_2, alpha_2, beta_2],
                       [beta_2, beta_2, alpha_2]]
            self.omega = [0.2250000000, 0.1323941527, 0.1323941527, 0.1323941527,
                          0.1259391805, 0.1259391805, 0.1259391805, 0.1259391805]

    def integrate_2d(self, func, meas, x1, x2, x3):
        f_int = 0
        # gauss
        for i in xrange(self.n):
            f_int += meas * self.omega[i] * \
                     func(self.xi[i][0] * x1[0] + self.xi[i][1] * x2[0] + self.xi[i][2] * x3[0],
                          self.xi[i][0] * x1[1] + self.xi[i][1] * x2[1] + self.xi[i][2] * x3[1])
        return f_int

class SpaceIntegration3d():
    def __init__(self, n):
        # order
        self.n = n
        # gauss
        if n == 1:
            self.omega = [0.25, 0.25, 0.25, 0.25]
            self.omega = [1.0]
        elif n == 4:
            alpha = 0.58541020
            beta = 0.13819660
            self.xi = [[alpha, beta, beta, beta],
                       [beta, alpha, beta, beta],
                       [beta, beta, alpha, beta],
                       [beta, beta, beta, alpha]]
            self.omega = [0.25, 0.25, 0.25, 0.25]
        elif n == 5:
            alpha = 0.5
            beta = 1.0 / 6.0
            self.xi = [[0.25, 0.25, 0.25, 0.25],
                       [alpha, beta, beta, beta],
                       [beta, alpha, beta, beta],
                       [beta, beta, alpha, beta],
                       [beta, beta, beta, alpha]]
            self.omega = [-16.0 / 20.0, 9.0 / 20.0, 9.0 / 20.0, 9.0 / 20.0, 9.0 / 20.0]

    def integrate_3d(self, func, volume, x1, x2, x3, x4):
        f_int = 0
        # gauss
        for i in xrange(self.n):
            f_int += volume * self.omega[i] * \
                     func(self.xi[i][0] * x1[0] + self.xi[i][1] * x2[0] + self.xi[i][2] * x3[0] + self.xi[i][3] * x4[0],
                          self.xi[i][0] * x1[1] + self.xi[i][1] * x2[1] + self.xi[i][2] * x3[1] + self.xi[i][3] * x4[1],
                          self.xi[i][0] * x1[2] + self.xi[i][1] * x2[2] + self.xi[i][2] * x3[2] + self.xi[i][3] * x4[2])
        return f_int
