__author__ = 'svetlana'

from math import sqrt

class Space1dIntegration():
    def __init__(self, n):
        self.n = n
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

class Space2dIntegration():
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


class SpaceIntegrator():
    def __init__(self, n, dim):
        # problem dimension
        self.dim = dim
        # order
        self.n = n

        # Depending on dimension of the problem, define points and weights
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
                f_int += 0.5 * measure * self.omega[i] * \
                         func(self.xi[i][0] * x_n[0, 0] + self.xi[i][1] * x_n[0, 1] + self.xi[i][2] * x_n[0, 2],
                              self.xi[i][0] * x_n[1, 0] + self.xi[i][1] * x_n[1, 1] + self.xi[i][2] * x_n[1, 2])
            elif self.dim == 3:
                f_int += measure * self.omega[i] * \
                         func(self.xi[i][0] * x_n[0, 0] + self.xi[i][1] * x_n[0, 1] + self.xi[i][2] * x_n[0, 2] + self.xi[i][3] * x_n[0, 3],
                              self.xi[i][0] * x_n[1, 0] + self.xi[i][1] * x_n[1, 1] + self.xi[i][2] * x_n[1, 2] + self.xi[i][3] * x_n[1, 3],
                              self.xi[i][0] * x_n[2, 0] + self.xi[i][1] * x_n[2, 1] + self.xi[i][2] * x_n[2, 2] + self.xi[i][3] * x_n[2, 3])
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
            self.omega = [-4.0 / 5.0, 9.0 / 20.0, 9.0 / 20.0, 9.0 / 20.0, 9.0 / 20.0]

    def integrate_3d(self, func, volume, x1, x2, x3, x4):
        f_int = 0
        # gauss
        for i in xrange(self.n):
            f_int += volume * self.omega[i] * \
                     func(self.xi[i][0] * x1[0] + self.xi[i][1] * x2[0] + self.xi[i][2] * x3[0] + self.xi[i][3] * x4[0],
                          self.xi[i][0] * x1[1] + self.xi[i][1] * x2[1] + self.xi[i][2] * x3[1] + self.xi[i][3] * x4[1],
                          self.xi[i][0] * x1[2] + self.xi[i][1] * x2[2] + self.xi[i][2] * x3[2] + self.xi[i][3] * x4[2])
        return f_int
