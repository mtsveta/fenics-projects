__author__ = 'svetlana'

from dolfin import *

class TimeIntegrator():
    def __init__(self, n, t_k, t_k1):
        self.tau = 0
        self.v = [0.0 for x in range(n)]
        self.n = n
        # gauss
        if n == 2:
            self.xi = [-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)]
            self.omega = [1.0, 1.0]
        elif n == 3:
            self.xi = [-sqrt(0.6), 0, sqrt(0.6)]
            self.omega = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
        elif n == 4:
            self.xi = [-0.861136311594953, -0.339981043584856, 0.339981043584856, 0.861136311594953]
            self.omega = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]

    def update(self, t_k, t_k1):
        self.tau = t_k1 - t_k
        # gauss
        for i in xrange(self.n):
            self.v[i] = 0.5 * (t_k * (1.0 - self.xi[i]) + t_k1 * (1.0 + self.xi[i]))

    __update = update  # private copy of original update() method

    def f(self, func, V):
        f_int = 0
        for i in xrange(self.n):
            func.t = self.v[i]
            c_i = 0.5 * self.tau * self.omega[i]
            f_int += interpolate(func, V) * c_i

        return f_int

    def sum_functions(self, f_vals):
        if self.n == 4:
            f_int = f_vals[0] + f_vals[1] + f_vals[2] + f_vals[3]
        elif self.n == 3:
            f_int = f_vals[0] + f_vals[1] + f_vals[2]
        elif self.n == 2:
            f_int = f_vals[0] + f_vals[1]

        return f_int

    def f_t_tk_tk1_t(self, func, V):
        f_int = 0
        for i in xrange(self.n):
            func.t = self.v[i]
            c_i = 0.125 * (self.tau ** 3) * (1.0 - self.xi[i] ** 2) * self.omega[i]
            f_int += interpolate(func, V) * c_i
        return f_int

    def f_t_tk(self, func, V):
        f_vals = []
        for i in xrange(self.n):
            func.t = self.v[i]
            c_i = 0.25 * (self.tau ** 2) * (1.0 + self.xi[i]) * self.omega[i]
            f_vals.append(interpolate(func, V) * c_i)

        f_int = self.sum_functions(f_vals)
        '''
        if self.n == 4:
            f_int = f_vals[0] + f_vals[1] + f_vals[2] + f_vals[3]
        elif self.n == 3:
            f_int = f_vals[0] + f_vals[1] + f_vals[2]
        elif self.n == 2:
            f_int = f_vals[0] + f_vals[1]
        '''

        return f_int

    def f_tk1_t(self, func, V):
        f_vals = []
        for i in xrange(self.n):
            func.t = self.v[i]
            c_i = 0.25 * (self.tau ** 2) * (1.0 - self.xi[i]) * self.omega[i]
            f_vals.append(interpolate(func, V) * c_i)

        f_int = self.sum_functions(f_vals)
        '''
        if self.n == 4:
            f_int = f_vals[0] + f_vals[1] + f_vals[2] + f_vals[3]
        elif self.n == 3:
            f_int = f_vals[0] + f_vals[1] + f_vals[2]
        elif self.n == 2:
            f_int = f_vals[0] + f_vals[1]
        '''
        return f_int

class SpaceIntegrator():
    def __init__(self, n, dim):
        # problem dimension
        self.dim = dim
        # order
        self.n = n
        if dim == 1:
            # gauss weight for coordination points xi and weights for function in this points
            if n == 2:
                self.xi = [-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)]
                self.omega = [1.0, 1.0]
            elif n == 3:
                self.xi = [-sqrt(0.6), 0, sqrt(0.6)]
                self.omega = [5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0]
            elif n == 4:
                self.xi = [-0.861136311594953, -0.339981043584856, 0.339981043584856, 0.861136311594953]
                self.omega = [0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454]
        elif dim == 2:
            # gauss weight for coordination points xi and weights for function in this points
            if n == 3:
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
        elif dim == 3:
            # gauss weight for coordination points xi and weights for function in this points
            if n == 4:
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

    def integrate(self, func, measure, x_n):
        f_int = 0
        # gauss
        for i in xrange(self.n):
            # x_{dim, dim + 1}
            if self.dim == 1:
                f_int += 0.5 * measure * self.omega[i] * \
                         func(0.5 * (x_n[0, 0] * (1 - self.xi[i]) + x_n[0, 1] * (1 + self.xi[i])))
            elif self.dim == 2:
                f_int += measure * self.omega[i] * \
                         func(self.xi[i][0] * x_n[0, 0] + self.xi[i][1] * x_n[0, 1] + self.xi[i][2] * x_n[0, 2],
                              self.xi[i][0] * x_n[1, 0] + self.xi[i][1] * x_n[1, 1] + self.xi[i][2] * x_n[1, 2])
            elif self.dim == 3:
                f_int += measure * self.omega[i] * \
                         func(self.xi[i][0] * x_n[0, 0] + self.xi[i][1] * x_n[0, 1] + self.xi[i][2] * x_n[0, 2] + self.xi[i][3] * x_n[0, 3],
                              self.xi[i][0] * x_n[1, 0] + self.xi[i][1] * x_n[1, 1] + self.xi[i][2] * x_n[1, 2] + self.xi[i][3] * x_n[1, 3],
                              self.xi[i][0] * x_n[2, 0] + self.xi[i][1] * x_n[2, 1] + self.xi[i][2] * x_n[2, 2] + self.xi[i][3] * x_n[2, 3])
        return f_int
