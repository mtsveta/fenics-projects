__author__ = 'svetlana'

from dolfin.cpp.common import tic, toc
from dolfin import *
import numpy
import math
from dolfin.cpp.mesh import CellFunction, cells
import integrators
import postprocess

def Div(f, dim):
    # Define differential operators for the 1D case
    if dim > 1:
        div_f = div(f)
    else:
        div_f = Dx(f, 0)
    return div_f

def get_matrices_of_optimization_problem(H_div, F_t_tk_, v_k, v_k1, grad_v_k, grad_v_k1, mesh, tau, dim):
    # Define variational problem to optimize the majorant with respect to y_1
    y = TrialFunction(H_div)
    q = TestFunction(H_div)

    S = assemble(Div(y, dim) * Div(q, dim) * dx(domain=mesh))
    K = assemble(inner(y, q) * dx(domain=mesh))

    # if w(x, t) = w(x) * (t - t_k) / tau
    z = assemble((F_t_tk_ -  Constant(0.5) * tau * (v_k1 - v_k)) * Div(q, dim) * dx(domain=mesh))
    q = assemble(inner(Constant(0.5) * grad_v_k + grad_v_k1, q) * dx(domain=mesh))

    return S, K, z, q

def update_majorant_on_k1_time_level(mf_y_k1, md_y_k1):
    Mf_y_k1 = assemble(mf_y_k1)
    Md_y_k1 = assemble(md_y_k1)
    return Mf_y_k1, Md_y_k1

def update_majorant_components(Md_y_k, Md_y_k1, Mf_y_k, Mf_y_k1):
    m_d = Md_y_k + Md_y_k1
    m_f = Mf_y_k + Mf_y_k1
    return m_d, m_f

def output_optimization_results(iter, maj, m_d, m_f, i_eff, beta, e_incr_k, ed_k, et_k):
    print 'opt cycle # %d' % iter
    print 'maj   = %8.2e, m_d = %8.2e, m_f = %8.2e, i_eff = %.4f, beta = %.4f' \
          % (maj, m_d, m_f, i_eff, beta)
    print 'error = %8.2e, e_d = %8.2e, e_T = %8.2e' % (e_incr_k, ed_k, et_k)

def calculate_majorant(m_d, m_f, C_FD, delta):
    # Calculate the increment of majorant on the time-cylinder [t_k, t_{k+1}]
    if m_f <= DOLFIN_EPS:
        beta = 1.0
        maj = 1.0 / delta * m_d
    else:
       if m_d <= DOLFIN_EPS:
           beta = 1.0
           maj = 1.0 / delta * C_FD ** 2 * m_f
       else:
           beta = C_FD * sqrt(m_f / m_d)
           maj = 1.0 / delta * ((1.0 + beta) * m_d + (1.0 + 1.0 / beta) * C_FD ** 2 * m_f)

    return maj, beta

def majorant_nd_t(k, m_incr_k, m0_k, md_k, mf_k, beta_k,
                  e_incr_k, ed_k, et_k,
                  f, sq_f, phi, quadr,
                  v_k, v_k1, grad_v_k, grad_v_k1, V_exact,
                  y_k, y_k1, H_div,
                  mesh, tau, delta, dim, domain, C_FD,
                  majorant_reconstruction_time,
                  majorant_minimization_time):
    tic()

    rf_k = Div(y_k, dim) - (v_k1 - v_k) / tau
    rf_k1 = Div(y_k1, dim) - (v_k1 - v_k) / tau

    rd_k = y_k - grad_v_k
    rd_k1 = y_k1 - grad_v_k1

    sq_F_ = quadr.f(sq_f, V_exact)
    F_t_tk_ = quadr.f_t_tk(f, V_exact)
    F_tk1_t_ = quadr.f_tk1_t(f, V_exact)

    # Components of mf = || f + div y - v_t ||
    mf_y_k = (sq_F_ + Constant(2.0) / tau * inner(F_tk1_t_, rf_k) + tau / Constant(3.0) * inner(rf_k, rf_k)) * dx(domain=mesh)
    mf_y_k1 = (Constant(2.0) / tau * inner(F_t_tk_, rf_k1) + tau / Constant(3.0) * (inner(rf_k, rf_k1) + inner(rf_k1, rf_k1))) * dx(domain=mesh)

    # Components of md = || y - grad v ||
    md_y_k = (tau / Constant(3.0) * inner(rd_k, rd_k)) * dx(domain=mesh)
    md_y_k1 = (tau / Constant(3.0) * (inner(rd_k, rd_k1) + inner(rd_k1, rd_k1))) * dx(domain=mesh)

    m0 = assemble(((phi - v_k) * (phi - v_k)) * dx(domain=mesh))

    # Assemble components of dual and balance terms of majorant
    # which stays the same in optimization process
    Mf_y_k = assemble(mf_y_k)
    Md_y_k = assemble(md_y_k)

    Mf_y_k1, Md_y_k1 = update_majorant_on_k1_time_level(mf_y_k1, md_y_k1)
    m_d, m_f = update_majorant_components(Md_y_k, Md_y_k1, Mf_y_k, Mf_y_k1)

    maj, beta = calculate_majorant(m_d, m_f, C_FD, delta)
    i_eff = sqrt(maj / e_incr_k)

    majorant_reconstruction_time[k] = toc()
    output_optimization_results(-1, maj, m_d, m_f, i_eff, beta, e_incr_k, ed_k, et_k)

    #----------------------------------------------------------------------------#
    # Optimization part
    #----------------------------------------------------------------------------#
    tic()
    S, K, z, q = get_matrices_of_optimization_problem(H_div, F_t_tk_, v_k, v_k1, grad_v_k, grad_v_k1, mesh, tau, dim)

    OPT_ITER = 5
    for opt_iter in range(1, OPT_ITER):
        # Solve system with respect to Y
        A = (C_FD ** 2) * S + beta * K
        b = - 0.5 * A * y_k.vector() - C_FD ** 2 * 3.0 / (tau ** 2) * z + beta * q
        #solve(A, y_k1.vector(), b, "gmres", 'ilu')
        solve(A, y_k1.vector(), b)

        Mf_y_k1, Md_y_k1 = update_majorant_on_k1_time_level(mf_y_k1, md_y_k1)
        m_d, m_f = update_majorant_components(Md_y_k, Md_y_k1, Mf_y_k, Mf_y_k1)

        maj, beta = calculate_majorant(m_d, m_f, C_FD, delta)
        i_eff = sqrt(maj / e_incr_k)
        output_optimization_results(opt_iter, maj, m_d, m_f, i_eff, beta, e_incr_k, ed_k, et_k)
    majorant_minimization_time[k] = toc()

    md_k[k] = m_d
    mf_k[k] = m_f
    beta_k[k] = beta
    m_incr_k[k] = maj
    m0_k[k] = m0

    return m_incr_k, m0_k, md_k, mf_k, beta_k, y_k1

def allocate_array_1d(n):
    return numpy.array([0.0 for i in range(n)])

def allocate_array_2d(m, n):
    return numpy.array([[0.0 for j in range(n)] for i in range(m)])

def error_majorant_distribution_nd(mesh, dim,
                                   grad_u_e, sq_grad_u_e, quad,
                                   grad_v_k, grad_v_k1, V_exact, VV_exact,
                                   y_k, y_k1, tau):
    cell_num = mesh.num_cells()
    e_distr = allocate_array_1d(cell_num)
    m_distr = allocate_array_1d(cell_num)

    rd_k = y_k - grad_v_k
    rd_k1 = y_k1 - grad_v_k1

    sq_grad_u = quad.f(sq_grad_u_e, V_exact)
    grad_u_t_t0 = quad.f_t_tk(grad_u_e, VV_exact)
    grad_u_t1_t = quad.f_tk1_t(grad_u_e, VV_exact)

    md = project(tau / 3.0 * (inner(rd_k, rd_k) + inner(rd_k, rd_k1) + inner(rd_k1, rd_k1)), V_exact)
    ed = project(sq_grad_u
                 - 2.0 / tau * (inner(grad_u_t1_t, grad_v_k) + inner(grad_u_t_t0, grad_v_k1))
                 + tau / 3.0 * (inner(grad_v_k, grad_v_k)
                                + inner(grad_v_k1, grad_v_k)
                                + inner(grad_v_k1, grad_v_k1)), V_exact)

    dofmap = V_exact.dofmap()

    '''
    scheme_order = 4
    gauss = integrators.SpaceIntegrator(scheme_order, dim)
    for c in cells(mesh):
        # Obtaining the coordinates of the vertices in the cell
        verts = c.get_vertex_coordinates().reshape(dim+1, dim)
        x_n_transpose = verts.T

        # Calculating the area of the cell
        matrix = postprocess.allocate_array_2d(dim+1, dim+1)
        matrix[0:dim, 0:dim + 1] = x_n_transpose
        matrix[dim, 0:dim+1] = numpy.array([1.0 for i in range(dim + 1)])
        #print "matrix = ", matrix
        meas = abs(numpy.linalg.det(matrix))


        e_distr[c.index()] = gauss.integrate(ed, meas, x_n_transpose)
        m_distr[c.index()] = gauss.integrate(md, meas, x_n_transpose)
    '''

    if dim == 1:
        gauss = SpaceIntegration1d(4)
        for c in cells(mesh):
            verts = dofmap.tabulate_coordinates(c)
            x_1 = verts[0, 0]
            x_2 = verts[1, 0]
            e_distr[c.index()] = gauss.integrate_1d(ed, c.diameter(), x_1, x_2)
            m_distr[c.index()] = gauss.integrate_1d(md, c.diameter(), x_1, x_2)

    elif dim == 2:

        gauss = SpaceIntegration2d(4)
        for c in cells(mesh):
            verts = dofmap.tabulate_coordinates(c)
            x_1 = verts[0, :]
            x_2 = verts[1, :]
            x_3 = verts[2, :]
            matrix = [[x_1[0], x_2[0], x_3[0]],
                      [x_1[1], x_2[1], x_3[1]],
                      [1, 1, 1]]
            meas = 0.5 * abs(numpy.linalg.det(matrix))

            e_distr[c.index()] = gauss.integrate_2d(ed, meas, x_1, x_2, x_3)
            m_distr[c.index()] = gauss.integrate_2d(md, meas, x_1, x_2, x_3)

    elif dim == 3:
        gauss = SpaceIntegration3d(5)

        for c in cells(mesh):
            verts = dofmap.tabulate_coordinates(c)
            x_1 = verts[0, :]
            x_2 = verts[1, :]
            x_3 = verts[2, :]
            x_4 = verts[3, :]
            matrix = [[x_1[0], x_2[0], x_3[0], x_4[0]],
                      [x_1[1], x_2[1], x_3[1], x_4[1]],
                      [x_1[2], x_2[2], x_3[2], x_4[2]],
                      [1.0, 1.0, 1.0, 1.0]]
            meas = abs(numpy.linalg.det(matrix))

            e_distr[c.index()] = gauss.integrate_3d(ed, meas, x_1, x_2, x_3, x_4)
            m_distr[c.index()] = gauss.integrate_3d(md, meas, x_1, x_2, x_3, x_4)

    return e_distr, m_distr, ed, md

class SpaceIntegration1d():

    def __init__(self, n):
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

def add_increment_of_error_and_majorant(k, e_incr_k, maj_incr_k, error_k, et_k, majorant_k, i_eff_maj_k, v_norm_incr_k, v_norm_k, vt_k, rel_error_k, rel_majorant_k):
    if k == 0:  # first time interval
        v_norm_k[k] = v_norm_incr_k[k]
        error_k[k] = e_incr_k[k]
        majorant_k[k] = maj_incr_k[k]
    else:
        v_norm_k[k] = v_norm_k[k - 1] - vt_k[k - 1] + v_norm_incr_k[k]
        error_k[k] = error_k[k - 1] - et_k[k - 1] + e_incr_k[k]
        majorant_k[k] = majorant_k[k - 1] + maj_incr_k[k]


    rel_error_k[k] = error_k[k] / v_norm_k[k]
    rel_majorant_k[k] = majorant_k[k] / v_norm_k[k]
    i_eff_maj_k[k] = sqrt(majorant_k[k] / error_k[k])

    return error_k, majorant_k, i_eff_maj_k, rel_error_k, rel_majorant_k


def bulk_marking(mesh, distr, theta):
    marking = CellFunction("bool", mesh)
    distr_0 = sorted(distr, reverse=True)[int(len(distr) * theta)]
    for cell in cells(mesh):
        marking[cell] = distr[cell.index()] > distr_0

    return marking


def predefined_amount_of_elements_marking(mesh, distr, theta):

    cells_num = mesh.num_cells()
    marking = CellFunction("bool", mesh)
    marking.set_all(False)

    i_cut = int(math.floor(theta * cells_num))
    index_sorted = sorted(range(len(distr)), key=lambda k: distr[k], reverse=True)
    cells_indices_to_refine = index_sorted[0:i_cut]

    for cell in cells(mesh):
        if cell.index() in cells_indices_to_refine:
            marking[cell] = True

    return marking


def averaged_marking(mesh, distr):
    cells_num = mesh.num_cells()
    marking = CellFunction("bool", mesh)
    distr_aver = sum(distr) / cells_num
    for c in cells(mesh):
        marking[c] = distr[c.index()] >= distr_aver
        #print distr[c.index()] >= distr_aver
    return marking




def allocate_space_for_error_and_majorant(n_t):
    ed_k = numpy.array([0.0 for x in range(n_t)])
    et_k = numpy.array([0.0 for x in range(n_t)])
    e_incr_k = numpy.array([0.0 for x in range(n_t)])
    error_k = numpy.array([0.0 for x in range(n_t)])

    md_k = numpy.array([0.0 for x in range(n_t)])
    mf_k = numpy.array([0.0 for x in range(n_t)])
    m_incr_k = numpy.array([0.0 for x in range(n_t)])
    majorant_k = numpy.array([0.0 for x in range(n_t)])

    beta_k = numpy.array([0.0 for x in range(n_t)])
    i_eff_mak_k = numpy.array([0.0 for x in range(n_t)])

    return ed_k, et_k, e_incr_k, e_incr_k, error_k, md_k, mf_k, m_incr_k, majorant_k, beta_k, i_eff_mak_k

def majorant_II_nd_t(k, m_incr_k, m0_k, md_k, mf_k, beta_k,
                     e_incr_k, ed_k, et_k,
                     f, sq_f, phi, quadr,
                     v_k, v_k1, grad_v_k, grad_v_k1, V_exact,
                     y_k, y_k1, H_div,
                     mesh, tau, delta, dim, domain, C_FD,
                     majorant_reconstruction_time,
                     majorant_minimization_time):
    tic()

    rf_k = Div(y_k, dim) - (v_k1 - v_k) / tau
    rf_k1 = Div(y_k1, dim) - (v_k1 - v_k) / tau

    rd_k = y_k - grad_v_k
    rd_k1 = y_k1 - grad_v_k1

    sq_F_ = quadr.f(sq_f, V_exact)
    F_t_tk_ = quadr.f_t_tk(f, V_exact)
    F_tk1_t_ = quadr.f_tk1_t(f, V_exact)

    # Components of mf = || f + div y - v_t ||
    mf_y_k = (sq_F_ + Constant(2.0) / tau * inner(F_tk1_t_, rf_k) + tau / Constant(3.0) * inner(rf_k, rf_k)) * dx(domain=mesh)
    mf_y_k1 = (Constant(2.0) / tau * inner(F_t_tk_, rf_k1) + tau / Constant(3.0) * (inner(rf_k, rf_k1) + inner(rf_k1, rf_k1))) * dx(domain=mesh)

    # Components of md = || y - grad v ||
    md_y_k = (tau / Constant(3.0) * inner(rd_k, rd_k)) * dx(domain=mesh)
    md_y_k1 = (tau / Constant(3.0) * (inner(rd_k, rd_k1) + inner(rd_k1, rd_k1))) * dx(domain=mesh)

    m0 = assemble(((phi - v_k) * (phi - v_k)) * dx(domain=mesh))

    # Assemble components of dual and balance terms of majorant
    # which stays the same in optimization process
    Mf_y_k = assemble(mf_y_k)
    Md_y_k = assemble(md_y_k)

    Mf_y_k1, Md_y_k1 = update_majorant_on_k1_time_level(mf_y_k1, md_y_k1)
    m_d, m_f = update_majorant_components(Md_y_k, Md_y_k1, Mf_y_k, Mf_y_k1)

    maj, beta = calculate_majorant(m_d, m_f, C_FD, delta)
    i_eff = sqrt(maj / e_incr_k)

    majorant_reconstruction_time[k] = toc()
    output_optimization_results(-1, maj, m_d, m_f, i_eff, beta, e_incr_k, ed_k, et_k)

    #----------------------------------------------------------------------------#
    # Optimization part
    #----------------------------------------------------------------------------#
    tic()
    S, K, z, q = get_matrices_of_optimization_problem(H_div, F_t_tk_, v_k, v_k1, grad_v_k, grad_v_k1, mesh, tau, dim)

    OPT_ITER = 5
    for opt_iter in range(1, OPT_ITER):
        # Solve system with respect to Y
        A = (C_FD ** 2) * S + beta * K
        b = - 0.5 * A * y_k.vector() - C_FD ** 2 * 3.0 / (tau ** 2) * z + beta * q
        #solve(A, y_k1.vector(), b, "gmres", 'ilu')
        solve(A, y_k1.vector(), b)

        Mf_y_k1, Md_y_k1 = update_majorant_on_k1_time_level(mf_y_k1, md_y_k1)
        m_d, m_f = update_majorant_components(Md_y_k, Md_y_k1, Mf_y_k, Mf_y_k1)

        maj, beta = calculate_majorant(m_d, m_f, C_FD, delta)
        i_eff = sqrt(maj / e_incr_k)
        output_optimization_results(opt_iter, maj, m_d, m_f, i_eff, beta, e_incr_k, ed_k, et_k)
    majorant_minimization_time[k] = toc()

    md_k[k] = m_d
    mf_k[k] = m_f
    beta_k[k] = beta
    m_incr_k[k] = maj
    m0_k[k] = m0

    return m_incr_k, m0_k, md_k, mf_k, beta_k, y_k1
