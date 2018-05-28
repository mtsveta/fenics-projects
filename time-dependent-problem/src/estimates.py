__author__ = 'svetlana'

from dolfin.cpp.common import tic, toc
from dolfin import *
import numpy
import math
from dolfin.cpp.mesh import CellFunction, cells
import integrators
import postprocess
import problem


def get_matrices_of_optimization_problem(H_div, F_t_tk_, v_k, v_k1, grad_v_k, grad_v_k1, mesh, tau, dim):
    # Define variational problem to optimize the majorant with respect to y_1
    y = TrialFunction(H_div)
    q = TestFunction(H_div)

    S = assemble(problem.Div(y, dim) * problem.Div(q, dim) * dx(domain=mesh))
    K = assemble(inner(y, q) * dx(domain=mesh))

    # if w(x, t) = w(x) * (t - t_k) / tau
    z = assemble((F_t_tk_ - Constant(0.5) * tau * (v_k1 - v_k)) * problem.Div(q, dim) * dx(domain=mesh))
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

    rf_k = problem.Div(y_k, dim) - (v_k1 - v_k) / tau
    rf_k1 = problem.Div(y_k1, dim) - (v_k1 - v_k) / tau

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

    OPT_ITER = 4
    for opt_iter in range(1, OPT_ITER):

        # Solve system with respect to Y
        A = (C_FD ** 2) * S + beta * K
        b = - 0.5 * A * y_k.vector() - C_FD ** 2 * 3.0 / (tau ** 2) * z + beta * q
        #solve(A, y_k1.vector(), b, "gmres", 'ilu')
        solve(A, y_k1.vector(), b)

        mf_y_k1 = (Constant(2.0) / tau * inner(F_t_tk_, rf_k1)
                   + tau / Constant(3.0) * (inner(rf_k, rf_k1) + inner(rf_k1, rf_k1))) * dx(domain=mesh)
        md_y_k1 = (tau / Constant(3.0) * (inner(rd_k, rd_k1) + inner(rd_k1, rd_k1))) * dx(domain=mesh)
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


    rf_k = problem.Div(y_k, dim) - (v_k1 - v_k) / tau
    rf_k1 = problem.Div(y_k1, dim) - (v_k1 - v_k) / tau

    rd_k = y_k - grad_v_k
    rd_k1 = y_k1 - grad_v_k1

    md_var = tau / Constant(3.0) * (inner(rd_k, rd_k) + inner(rd_k, rd_k1) + inner(rd_k1, rd_k1))
    mf_var = 1.0 / delta * ((1.0 + beta) * md_var +
                            (1.0 + 1.0 / beta) * (C_FD ** 2) *
                            (sq_F_ + Constant(2.0) / tau * (inner(F_tk1_t_, rf_k) + inner(F_t_tk_, rf_k1))
                             + tau / Constant(3.0) * (inner(rf_k, rf_k) + inner(rf_k, rf_k1) + inner(rf_k1, rf_k1))))

    return m_incr_k, m0_k, md_k, mf_k, beta_k, y_k1, md_var, mf_var

def allocate_array_1d(n):
    return numpy.array([0.0 for i in range(n)])

def allocate_array_2d(m, n):
    return numpy.array([[0.0 for j in range(n)] for i in range(m)])


def majorant_distribution_DG0(mesh, ed_var, md_var, mdf_var):

    # Define the functional space used for distribution over cells
    DG0 = FunctionSpace(mesh, "DG", 0)
    w = TestFunction(DG0)

    e_d_DG0 = Function(DG0)
    m_d_DG0 = Function(DG0)
    maj_DG0 = Function(DG0)

    e_d_var = w * ed_var * dx(domain=mesh)
    m_d_var = w * md_var * dx(domain=mesh)
    m_df_var = w * mdf_var * dx(domain=mesh)

    # Assemble the variation form and dumping obtained vector into function from DG0
    assemble(e_d_var, tensor=e_d_DG0.vector())
    assemble(m_d_var, tensor=m_d_DG0.vector())
    assemble(m_df_var, tensor=maj_DG0.vector())

    ed_distr = e_d_DG0.vector().array()
    md_distr = m_d_DG0.vector().array()
    maj_distr = maj_DG0.vector().array()

    print "ed_distr  DG0 total = %8.2e" % numpy.sum(ed_distr)
    print "md_distr  DG0 total = %8.2e" % numpy.sum(md_distr)
    print "maj_distr DG0 total = %8.2e\n" % numpy.sum(maj_distr)

    return ed_distr, md_distr, maj_distr

def error_majorant_distribution_nd(mesh, dim, ed_var, md_var, maj_var, V_exact):

    cell_num = mesh.num_cells()

    ed_distr = allocate_array_1d(cell_num)
    md_distr = allocate_array_1d(cell_num)
    maj_distr = allocate_array_1d(cell_num)

    ed = project(ed_var, V_exact)
    md = project(md_var, V_exact)
    maj = project(maj_var, V_exact)

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
        meas = float(1.0 / math.factorial(dim)) * abs(numpy.linalg.det(matrix))

        ed_distr[c.index()] = gauss.integrate(ed, meas, x_n_transpose)
        md_distr[c.index()] = gauss.integrate(md, meas, x_n_transpose)
        maj_distr[c.index()] = gauss.integrate(maj, meas, x_n_transpose)

    print "ed_distr  int total = %8.2e" % numpy.sum(ed_distr)
    print "md_distr  int total = %8.2e" % numpy.sum(md_distr)
    print "maj_distr int total = %8.2e\n" % numpy.sum(maj_distr)

    return ed_distr, md_distr, maj_distr

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

    rf_k = problem.Div(y_k, dim) - (v_k1 - v_k) / tau
    rf_k1 = problem.Div(y_k1, dim) - (v_k1 - v_k) / tau

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

    OPT_ITER = 4
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

def grad_error_norm(grad_u_e, sq_grad_u_e, quadr, V_exact, VV_exact, grad_vk, grad_vk1, tau, mesh):

    sq_grad_u = quadr.f(sq_grad_u_e, V_exact)
    grad_u_t_tk = quadr.f_t_tk(grad_u_e, VV_exact)
    grad_u_tk1_t = quadr.f_tk1_t(grad_u_e, VV_exact)

    # this should be projected on Vector Function Space in 2d
    error_var = (sq_grad_u
                - Constant(2.0) / tau * (inner(grad_u_tk1_t, grad_vk)
                                         + inner(grad_u_t_tk, grad_vk1))
                + tau / Constant(3.0) * (inner(grad_vk, grad_vk)
                                         + inner(grad_vk1, grad_vk)
                                         + inner(grad_vk1, grad_vk1)))
    e_d = assemble(error_var * dx(domain=mesh))
    return e_d, error_var


def error_norm_nd_t(k, e_incr_k, ed_k, et_k,
                    grad_u_e, sq_grad_u_e,
                    quadr,
                    u_k1, v_k1, grad_v_k, grad_v_k1, V_exact, VV_exact,
                    mesh, dim, tau, delta):


    v_k1_Ve = interpolate(v_k1, V_exact)

    grad_v_k_Ve = grad_v_k
    grad_v_k1_Ve = grad_v_k1
    # grad_v_k_Ve = project(grad_v_k, VV_exact)
    # grad_v_k1_Ve = project(grad_v_k1, VV_exact)
    #grad_v_k_Ve = interpolate(grad_v_k, VV_exact)
    #grad_v_k1_Ve = interpolate(grad_v_k1, VV_exact)


    e_d, ed_var = grad_error_norm(grad_u_e, sq_grad_u_e, quadr, V_exact, VV_exact, grad_v_k_Ve, grad_v_k1_Ve, tau, mesh)
    #e_d = my_grad_error_norm(grad_u_e, sq_grad_u_e, quadr, V_exact, VV_exact, grad_v_k, grad_v_k1, tau, mesh, dim)
    e_t = assemble((inner((u_k1 - v_k1_Ve), (u_k1 - v_k1_Ve))) * dx(domain=mesh))

    ed_k[k] = e_d
    et_k[k] = e_t
    # Calculate the error increment
    e_incr_k[k] = (2.0 - delta) * e_d + e_t

    print '%------------------------------------------------------------------------------------%'
    print '% Error '
    print '%------------------------------------------------------------------------------------%'
    print "\| grad e \|^2_Q = %8.2e" % e_d
    print "\| e \|^2_T      = %8.2e" % e_t
    #print "\| e \|^2_inf                    = %8.2e" % linf_e
    #print "\| (lmbd - 0.5 div a)^0.5 e \|^2 = %8.2e" % delta_e
    #print "\| lmbd^0.5 e \|^2               = %8.2e" % lambda_e
    #print "\| a e \|^2_{A^{-1}}             = %8.2e" % a_e

    return e_incr_k, ed_k, et_k, ed_var


def v_norm_nd_t(k, v_norm_incr_k, vd_k, vt_k, v_k, v_k1, V_exact, VV_exact, mesh, dim, tau, delta):

    v_k1_Ve = interpolate(v_k1, V_exact)
    v_k_Ve = interpolate(v_k, V_exact)

    #grad_v_k_Ve = project(grad_v_k, VV_exact)
    #grad_v_k1_Ve = project(grad_v_k1, VV_exact)
    #grad_v_k_Ve = interpolate(grad_v_k, VV_exact)
    #grad_v_k1_Ve = interpolate(grad_v_k1, VV_exact)
    #grad_v_k_Ve = grad_v_k
    #grad_v_k1_Ve = grad_v_k1


    vd_k[k] = assemble(tau / Constant(3.0) * (inner(problem.Grad(v_k_Ve, dim), problem.Grad(v_k_Ve, dim))
                                              + inner(problem.Grad(v_k_Ve, dim), problem.Grad(v_k1_Ve, dim))
                                              + inner(problem.Grad(v_k1_Ve, dim), problem.Grad(v_k1_Ve, dim)))* dx(domain=mesh))

    #vd_k[k] = assemble(tau / Constant(3.0) * (inner(grad_v_k_Ve, grad_v_k1_Ve)
    #                                          + inner(grad_v_k1_Ve, grad_v_k_Ve)
    #                                          + inner(grad_v_k1_Ve, grad_v_k1_Ve)) * dx(domain=mesh))
    vt_k[k] = assemble((inner(v_k1_Ve, v_k1_Ve)) * dx(domain=mesh))

    # Calculate the error increment
    v_norm_incr_k[k] = (2.0 - delta) * vd_k[k] + vt_k[k]

    return v_norm_incr_k, vd_k, vt_k
