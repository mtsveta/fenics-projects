from dolfin.cpp.common import tic, toc
from ufl.tensoralgebra import Inverse

__author__ = 'svetlana'

from dolfin import *
from dolfin.cpp.mesh import cells, CellFunction
import postprocess, integration
import numpy
import math
from problem import Grad, Div, inverse

def output_optimization_results(iter, maj, m_d, m_f, i_eff, beta, error):
    if iter >= 1:
        print 'opt cycle # %d' % iter
    print 'maj   = %8.2e, m_d = %8.2e, m_f = %8.2e, i_eff = %.4f, beta = %.4f' \
          % (maj, m_d, m_f, i_eff, beta)
    print 'error = %8.2e' % (error)

def calculate_majorant_bar_mu_opt(u, y, beta, f_bar, A, lambda_1, lmbd, mesh, dim, C_FD):
    """
    :param u: approximate solution
    :param y: flux
    :param beta: parameter minimizing majorant
    :param f: right-hand side function
    :param lmbd: reaction function
    :param mesh: mesh
    :param dim: problem dimension
    :param C_FD: Freidrichs constant
    :return maj: majorant value
    :return m_d, m_f_one_minus_mu_opt: majorant components
    :return beta: optimal parameter
    """
    # Define parameters
    mu_opt = (C_FD ** 2) * (1.0 + beta) * lmbd / (beta * lambda_1 + (C_FD ** 2) * (1.0 + beta) * lmbd)
    w_opt = (C_FD ** 2) * (1 + beta) / (beta * lambda_1 + (C_FD ** 2) * (1 + beta) * lmbd)

    # Define residuals
    r_d = y - A * Grad(u, dim)
    r_f = (Div(y, dim) + f_bar)

    # Define variational forms
    var_m_d = inner(inverse(A, dim) * r_d, r_d)
    var_m_f_w_opt = w_opt * inner(r_f, r_f)
    var_m_f_one_minus_mu_opt = ((1 - mu_opt) ** 2) * inner(r_f, r_f)

    # Define majorant components
    m_d = assemble(var_m_d * dx(domain=mesh))
    m_f_w_opt = assemble(var_m_f_w_opt * dx(domain=mesh))     # for calculating majorant
    m_f_one_minus_mu_opt = assemble(var_m_f_one_minus_mu_opt * dx(domain=mesh))    # fo calculating beta_opt

    # Calculate majorant based on the parameter value
    if m_f_w_opt <= DOLFIN_EPS:
        maj = m_d
    else:
       if m_d <= DOLFIN_EPS:
           maj = m_f_w_opt
       else:
           maj = (1.0 + beta) * m_d + m_f_w_opt

    # Calculate the optimal value for beta parameter
    beta = C_FD * sqrt(m_f_one_minus_mu_opt / m_d / lambda_1)

    return maj, m_d, m_f_one_minus_mu_opt, beta, var_m_d, var_m_f_w_opt

def get_matrices_of_optimization_problem_bar(H_div, v, f_bar, A, mesh, dim):

    # Define variational problem
    y = TrialFunction(H_div)
    q = TestFunction(H_div)

    # Define system of linear equation to find the majorant
    S = assemble(inner(Div(y, dim), Div(q, dim)) * dx(domain=mesh))
    K = assemble(inner(inverse(A, dim) * y, q) * dx(domain=mesh))
    z = assemble(inner(-f_bar, Div(q, dim)) * dx(domain=mesh))
    g = assemble(inner(Grad(v, dim), q) * dx(domain=mesh))

    return S, K, z, g


def majorant_nd(v, y, H_div, f, A, a, lambda_1, lmbd,
                error,
                mesh, dim, domain, C_FD,
                MAJORANT_OPTIMIZE):

    tic()

    # Initialize value
    beta = 1.0

    #for i in range(2):
    f_bar = f - lmbd * v - inner(a, Grad(v, dim))
    #f_bar = f
    maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = calculate_majorant_bar_mu_opt(v, y, beta, f_bar, A, lambda_1, lmbd, mesh, dim, C_FD)
    i_eff = sqrt(maj / error)
    output_optimization_results(-1, maj, m_d, m_f, i_eff, beta, error)

    majorant_reconstruction_time = toc()

    if MAJORANT_OPTIMIZE:
        #----------------------------------------------------------------------------#
        # Optimization algorithm
        #----------------------------------------------------------------------------#
        print " "
        print "%-----------------------"
        print "% Majorant optimization "
        print "%-----------------------"
        print " "

        tic()
        S, K, z, g = get_matrices_of_optimization_problem_bar(H_div, v, f_bar, A, mesh, dim)

        y = Function(H_div)
        Y = y.vector()

        OPT_ITER = 3
        # Execute iterative process to optimize majorant with respect to beta and flux
        for k in range(1, OPT_ITER):
            # Solve system with respect to Y
            solve((C_FD ** 2)/lambda_1 * S + beta * K, Y, (C_FD ** 2)/lambda_1 * z + beta * g)

            y.vector = Y
            # Calculate majorant
            maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = calculate_majorant_bar_mu_opt(v, y, beta, f_bar, A, lambda_1, lmbd, mesh, dim, C_FD)
            i_eff = sqrt(maj / error)

            output_optimization_results(k, maj, m_d, m_f, i_eff, beta, error)

        majorant_minimization_time = toc()

    else:

        majorant_minimization_time = 0.0

    return maj, y, beta, m_d, m_f, var_m_d, var_m_f_w_opt, majorant_reconstruction_time, majorant_minimization_time

def error_majorant_distribution_nd(mesh, dim,
                                   V_exact, var_grad_e, var_delta_e, var_m_d, var_m_f_w_opt, beta):

    # Allocate arrays for error and majorant distributions
    cell_num = mesh.num_cells()
    ed_distr = postprocess.allocate_array(cell_num)
    delta_e_distr = postprocess.allocate_array(cell_num)
    md_distr = postprocess.allocate_array(cell_num)
    mf_distr = postprocess.allocate_array(cell_num)

    # Project UFL forms on the high-order functional space to obtain functions
    ed = project(var_grad_e, V_exact)
    delta_e = project(var_delta_e, V_exact)
    md = project(var_m_d, V_exact)
    mf_wopt = project(var_m_f_w_opt, V_exact)

    # Obtained the mapping
    dofmap = V_exact.dofmap()

    # Define integration scheme order and assign the integrator
    scheme_order = 4
    gauss = integration.SpaceIntegrator(scheme_order, dim)

    # Iterate over cells of the mesh and obtain local value of the majorant
    for c in cells(mesh):

        # Obtaining the coordinates of the vertices in the cell
        verts = dofmap.tabulate_coordinates(c)
        x_n = verts[0:dim+1, :]
        x_n_transpose = x_n.transpose()

        # Calculating the area of the cell
        matrix = postprocess.allocate_array_2d(dim+1, dim+1)
        matrix[0:dim, 0:dim+1] = x_n_transpose
        matrix[dim, 0:dim+1] = numpy.array([1.0 for i in range(dim + 1)])
        meas = abs(numpy.linalg.det(matrix))

        # Integrating over the cell
        ed_distr[c.index()] = gauss.integrate(ed, meas, x_n_transpose)
        delta_e_distr[c.index()] = gauss.integrate(delta_e, meas, x_n_transpose)
        md_distr[c.index()] = gauss.integrate(md, meas, x_n_transpose)
        mf_distr[c.index()] = gauss.integrate(mf_wopt, meas, x_n_transpose)

    # Calculate majorant based on the parameter value
    if sum(mf_distr) <= DOLFIN_EPS:
        maj_distr = md_distr
    else:
       if sum(md_distr) <= DOLFIN_EPS:
           maj_distr = mf_distr
       else:
           maj_distr = (1.0 + beta) * md_distr + mf_distr

    e_distr = ed_distr + delta_e_distr

    '''
    print 'md_cells', md_distr
    print 'mf_cells', mf_distr
    print 'ed_cells', ed_distr
    print 'delta_e_cells', delta_e_distr

    print 'maj_cells', maj_distr
    print 'e_cells', e_distr

    print 'sum maj_cells', sum(maj_distr)
    print 'sum e_cells', sum(e_distr)
    '''

    return ed_distr, delta_e_distr, e_distr, md_distr, mf_distr, maj_distr, ed, md

def error_norm(u, ue, lmbd, A, a, Ve, mesh, dim):
    """
    :param u: approximate solution
    :param ue: exact solution
    :param lmbd: reaction function
    :param A: diffusion operator
    :param Ve: functional space of exact solution
    :param mesh: mesh
    :param dim: dimension of the problem

    :return: error-norm between u and ue
    """
    # Interpolate exact and approximate solution to the functional space of exact solution
    u_ve = interpolate(u, Ve)
    u_exact_ve = interpolate(ue, Ve)
    e = abs(u_ve - u_exact_ve)

    # Define variational form of the error
    var_grad_e = inner(A * Grad(e, dim), Grad(e, dim))
    delta = lmbd - 0.5 * Div(a, dim)
    var_delta_e = delta * inner(e, e)

    # Assembling variational form of the error
    grad_e = assemble(var_grad_e * dx(domain=mesh))
    delta_e = assemble(var_delta_e * dx(domain=mesh))

    print '%-----------------'
    print '% Error '
    print '%-----------------'
    print "\| grad e \|^2_A", grad_e
    print "\| (lmbd - 0.5 div a)^0.5 e \|^2", delta_e
    print "||| e |||^2", grad_e + delta_e

    return grad_e + delta_e, var_grad_e, var_delta_e

def minorant(v, delta, mesh, Vh, u0, u0_boundary, f, A, lmbd, a, dim, error):

    # Define variational problem
    w = TrialFunction(Vh)
    mu = TestFunction(Vh)

    # Define variational forms
    # a(w, mu)
    a = inner(A * Grad(w, dim), Grad(mu, dim) + lmbd * inner(w, mu) + inner(inner(a, Grad(w, dim)), mu)
              + delta * (inner(Div(Grad(w, dim), dim) + lmbd * w + inner(a, Grad(w, dim)), inner(a, Grad(mu, dim))))) * dx(domain=mesh)
    # F(u, mu)
    L = (f * mu - inner(A * Grad(v, dim), Grad(mu, dim)) - lmbd * inner(v, mu) - inner(inner(a, Grad(w, dim)), mu) +
         + delta * (inner(f, inner(a, Grad(mu, dim))))) * dx(domain=mesh)


    # Define boundary condition and solve the minimization problem with respect to w
    w = Function(Vh)
    bc = DirichletBC(Vh, Constant(0.0), u0_boundary)
    solve(a == L, w, bc)

    # Define variational form for minorant: 2 F(u, w) - a(w, w)
    min_var = (2 * (f * w - inner(A * Grad(v, dim), Grad(w, dim)) - lmbd * inner(v, w) - inner(inner(a, Grad(v, dim)), w))
               - inner(A * Grad(w, dim), Grad(w, dim)) - (lmbd - 0.5 * Div(a, dim)) * inner(w, w)) * dx(domain=mesh)
    min = assemble(min_var)

    return min, w

def output_result_error_and_majorant(error, majorant, i_eff_maj):
    print ' '
    print '%-----------------'
    print '% Result majorant '
    print '%-----------------'
    print 'error     = %8.4e ' % (error)
    print 'majorant  = %8.4e' % (majorant)
    print 'i_eff_maj = %.4f' % (i_eff_maj)
    print " "

def output_result_minorant(minorant, i_eff_min):
    print ' '
    print '%-----------------'
    print '% Result minorant '
    print '%-----------------'
    print 'minorant  = %8.4e' % (minorant)
    print 'i_eff_min = %.4f' % (i_eff_min)
    print " "

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

# goal functional - example 1 (green book)
def J_energy_norm_error(w, grad_e, norm_grad_e, dim):
    return inner(Grad(w, dim), grad_e) / norm_grad_e

def solve_dual_problem(mesh, V_star, u0, u0_boundary, u, u_e, dim, A, adjA, norm_grad_e):

    # Define dual variational problem
    z = TrialFunction(V_star)
    w = TestFunction(V_star)

    # Define the adjoint right-hand side
    a_star = inner(A * Grad(z, dim), Grad(w, dim)) * dx(domain=mesh)

    # Define the adjoint left-hand side
    L_star = J_energy_norm_error(w, A * Grad(u_e - u, dim), norm_grad_e, dim) * dx(domain=mesh)

    z = Function(V_star)
    bc_star = DirichletBC(V_star, u0, u0_boundary)
    solve(a_star == L_star, z, bc_star)

    z_exact = project( (u - u_e) / norm_grad_e, V_star)

    z_error = assemble((z - z_exact)**2 * dx(domain=mesh))
    print "z error", z_error

    return z


def get_indicators_CG0(mesh, V, V_star, f, A, adjA, lmbd, u0, u0_boundary, u, u_e, e_form, norm_grad_e, dim):

    # Contruct the solution of the adjoint problem
    z = solve_dual_problem(mesh, V_star, u0, u0_boundary, u, u_e, dim, A, adjA, norm_grad_e)
    Pz = project(z, V)

    # Get parameters of the mesh
    h = mesh.hmax()
    n = FacetNormal(mesh)

    # Construct the set of
    DG0 = FunctionSpace(mesh, "DG", 0)
    w = TestFunction(DG0)

    r_h = jump(A * Grad(u, dim), n)
    R_h = Div(A * Grad(u, dim), dim) + f - lmbd


    eta_var = w * (h * R_h)**2 * dx(domain=mesh) + avg(w) * avg(h) * r_h**2 * dS(domain=mesh)
    E_DWR_var = w * inner(R_h, z - Pz) * dx(domain=mesh) - 0.5 * avg(w) * inner(r_h, avg(z - Pz)) * dS(domain=mesh)
    J_e_var = w * J_energy_norm_error(e_form, Grad(e_form, dim), norm_grad_e, dim) * dx(domain=mesh)

    eta_DG0 = Function(DG0)
    E_DWR_DG0 = Function(DG0)
    J_e_DG0 = Function(DG0)

    assemble(eta_var, tensor=eta_DG0.vector())
    assemble(E_DWR_var, tensor=E_DWR_DG0.vector())
    assemble(J_e_var, tensor=J_e_DG0.vector())

    eta_distr = eta_DG0.vector().array()
    E_DWR_distr = E_DWR_DG0.vector().array()
    J_e_distr = J_e_DG0.vector().array()

    print "eta_DG0 total:", numpy.sum(eta_distr)
    print "E_DWR_DG0 total", numpy.sum(E_DWR_distr)
    print "J_e_DG0 total", numpy.sum(J_e_distr)

    return eta_distr, E_DWR_distr, J_e_distr

def majorant_distribution_DG0(mesh, f, lmbd, A, invA, u, e_form, y, beta, C_FD, dim):

    # Define the functional space used for distribution over cells
    DG0 = FunctionSpace(mesh, "DG", 0)
    w = TestFunction(DG0)

    m_d_DG0 = Function(DG0)
    m_df_DG0 = Function(DG0)
    e_d_DG0 = Function(DG0)

    # Define optimal parameters
    w_opt = (C_FD ** 2) * (1 + beta) / (beta + (C_FD ** 2) * (1 + beta) * lmbd)

    # Define the residuals of the majorant
    r_d = y - A * Grad(u, dim)
    r_f = (Div(y, dim) + f - lmbd * u)


    # Define variational forms of dual component of majorant and the whole functional
    m_d_var = w * sqrt(inner(invA * r_d, r_d)) * dx(domain=mesh)
    m_df_var = w * sqrt((1.0 + beta) * inner(invA * r_d, r_d) + w_opt * inner(r_f, r_f)) * dx(domain=mesh)
    e_d_var = w * sqrt(inner(A * Grad(e_form, dim), Grad(e_form, dim))) * dx(domain=mesh)

    # Assemble the variation form and dumping obtained vector into function from DG0
    assemble(m_d_var, tensor=m_d_DG0.vector())
    assemble(m_df_var, tensor=m_df_DG0.vector())
    assemble(e_d_var, tensor=e_d_DG0.vector())

    m_d_distr = m_d_DG0.vector().array()
    m_df_distr = m_df_DG0.vector().array()
    e_d_distr = e_d_DG0.vector().array()


    print "m_d_DG0 total", numpy.sum(m_d_distr)
    print "m_df_DG0 total", numpy.sum(m_df_distr)
    print "e_d_DG0 total", numpy.sum(e_d_distr)

    return m_d_distr, m_df_distr, e_d_distr

'''
# goal functional - example 1 (green book)
def J(w, grad_e, norm_grad_e, dim):
    return inner(Grad(w, dim), grad_e) / norm_grad_e

def solve_dual_problem(V_star, f, u0, u0_boundary, u, u_e, dim):

    # Define dual variational problem
    z = TrialFunction(V_star)
    w = TestFunction(V_star)

    # Define the system
    a_star = inner(Grad(z, dim), Grad(w, dim)) * dx(domain=mesh)
    grad_e = Grad(u_e, dim) - Grad(u, dim)
    norm_grad_e = sqrt(assemble(inner(grad_e, grad_e) * dx(domain = mesh)) )
    L_star = J(w, grad_e, norm_grad_e, dim)

    z = Function(V_star)
    bc_star = DirichletBC(V_star, u0, u0_boundary)
    solve(a_star == L_star, z, bc_star)

    z_exact = project( (u - u_e) / norm_grad_e, V_star)
    z_error = assemble((z - z_exact)**2 * dx(domain=mesh))

    print "z error", z_error

    return z

def compare_error_indicators(mesh, V, V_star, f, u0, u0_boundary, u, u_e, grad_u_e, y, beta, test_num, tag, dim):

    z = solve_dual_problem(V_star, f, u0, u0_boundary, u, u_e)
    Pz = project(z, V)

    norm_grad_e = sqrt(assemble(inner(grad_u_e - Grad(u, dim), grad_u_e - Grad(u, dim)) * dx(domain = mesh)))
    h = mesh.hmax()
    n = FacetNormal(mesh)
    DG0 = FunctionSpace(mesh, "DG", 0)
    w = TestFunction(DG0)

    r_h = jump(Grad(u, dim), n)
    R_h = div(Grad(u, dim)) + f

    cell_num = mesh.num_cells()

    e = abs(u_e - u)
    r_d = abs(Grad(u, dim) - y)
    r_f = abs(f + Div(y, dim))
    eta_var = w * (h * R_h)**2 * dx + avg(w) * avg(h) * r_h**2 * dS
    E_DWR_var = w * inner(R_h, z - Pz) * dx - 0.5 * avg(w) * inner(r_h, avg(z - Pz)) * dS
    J_e_var = w * J(w, e, norm_grad_e) * dx
    m_d_var = w * sqrt(inner(r_d, r_d)) * dx
    #C_FD = 1 / (sqrt(3) * DOLFIN_PI)
    height = 1
    width = 2
    length = 2
    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height**2 + 1.0 / width**2 + 1.0 / length**2)
    m_df_var = w * sqrt((1 + beta) * inner(r_d, r_d) + C_FD**2 * (1 + 1 / beta) * inner(r_f, r_f)) * dx

    #M_e = w * (u - u_exact) * dx

    eta_DG0 = Function(DG0)
    E_DWR_DG0 = Function(DG0)
    J_e_DG0 = Function(DG0)
    m_d_DG0 = Function(DG0)
    m_df_DG0 = Function(DG0)

    assemble(eta_var, tensor=eta_DG0.vector())
    assemble(E_DWR_var, tensor=E_DWR_DG0.vector())
    assemble(J_e_var, tensor=J_e_DG0.vector())
    assemble(m_d_var, tensor=m_d_DG0.vector())
    assemble(m_df_var, tensor=m_df_DG0.vector())

    eta_distr = eta_DG0.vector().array()
    E_DWR_distr = E_DWR_DG0.vector().array()
    J_e_distr = J_e_DG0.vector().array()
    m_d_distr = m_d_DG0.vector().array()
    m_df_distr = m_df_DG0.vector().array()

    eta_DG0_total = numpy.sum(eta_distr)
    E_DWR_DG0_total = numpy.sum(E_DWR_distr)
    J_e_DG0_total = numpy.sum(J_e_distr)
    m_d_DG0_total = numpy.sum(m_d_distr)
    m_df_DG0_total = numpy.sum(m_df_distr)

    print "eta_DG0 total:", eta_DG0_total
    print "E_DWR_DG0 total", E_DWR_DG0_total
    print "J_e_DG0 total", J_e_DG0_total
    print "m_d_DG0 total", m_d_DG0_total
    print "m_df_DG0 total", m_df_DG0_total

# goal functional - example 1 (green book)
def J_energy_norm_error(w, grad_e, norm_grad_e, dim):
    return inner(Grad(w, dim), grad_e) / norm_grad_e

def solve_dual_problem(mesh, V_star, u0, u0_boundary, e, dim, norm_grad_e):

    # Define dual variational problem
    z = TrialFunction(V_star)
    w = TestFunction(V_star)

    # Define the adjoint right-hand side
    a_star = inner(Grad(z, dim), Grad(w, dim)) * dx(domain=mesh)

    # Define the adjoint left-hand side
    L_star = J_energy_norm_error(w, Grad(e, dim), norm_grad_e, dim) * dx(domain=mesh)

    z = Function(V_star)
    bc_star = DirichletBC(V_star, u0, u0_boundary)
    solve(a_star == L_star, z, bc_star)

    z_exact = project( e / norm_grad_e, V_star)

    z_error = assemble((z - z_exact)**2 * dx(domain=mesh))
    print "z error", z_error

    return z
'''