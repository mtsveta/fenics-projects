__author__ = "Svetlana Matculevich <svetlana.matculevich@ricam..fi>"
__date__ = "2013-10-11"
__copyright__ = "Copyright (C) 2013 Svetlana Matculevich"
__license__ = ""

from dolfin import *
from dolfin.cpp.mesh import cells, CellFunction
from dolfin.cpp.common import tic, toc
import postprocess, integration
import numpy
import math
from problem import Grad, Div, inverse

# Console output of error and majorant values
def output_optimization_results(iter, maj, m_d, m_f, i_eff, beta, error):
    print 'opt. # %d: e^2 = %8.2e, maj^2 = %8.2e, m^2_d = %8.2e, m^2_f = %8.2e, (i_eff)^2 = %.4f' \
          % (iter, error, maj, m_d, m_f, i_eff)

# Calculate optimal (wrt to big changes in the reaction fucntion) majorant
def calculate_majorant_bar_mu_opt(u, y, beta, f_bar, A, invA, lambda_1, lmbd, mesh, dim, C_FD):
    """
    :param u: approximate solution
    :param y: flux
    :param beta: parameter minimizing majorant
    :param f_bar: right-hand side function (modified RHS)
    :param A: diffusion matrix
    :param lambda_1: minimal eigenvalue of diffusion matrix
    :param lmbd: reaction function
    :param mesh: mesh
    :param dim: geometrical problem dimension
    :param C_FD: Freidrichs constant
    :return maj: majorant value
            m_d, m_f_one_minus_mu_opt: majorant components
            beta: optimal parameter
            var_m_d, var_m_f_w_opt: variational form of the majorant (to use them later in construction of error estimators)
    """
    # Define optimal parameters
    mu_opt = (C_FD ** 2) * (1.0 + beta) * lmbd / (beta * lambda_1 + (C_FD ** 2) * (1.0 + beta) * lmbd)
    w_opt = (C_FD ** 2) * (1 + beta) / (beta * lambda_1 + (C_FD ** 2) * (1 + beta) * lmbd)

    # Define residuals
    r_d = y - A * Grad(u, dim)
    r_f = (Div(y, dim) + f_bar)

    # Define variational forms
    var_m_d = inner(invA * r_d, r_d)
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

def get_matrices_of_optimization_problem_bar(H_div, v, f_bar, invA, mesh, dim):
    """
    :param H_div: funtional space for the flux function
    :param v: approximate solution
    :param f_bar: right-hand side function (modified RHS)
    :param A: diffusion matrix
    :param mesh: mesh
    :param dim: geometrical problem dimension
    :return DivDiv, PhiPhi, RhsDivPhi, RhsPhi: matrixes for the optimal reconstruction of the flux
    """
    # Define variational problem
    y = TrialFunction(H_div)
    q = TestFunction(H_div)

    # Define system of linear equation to find the majorant
    # Phi_i, i = 1, ..., d are vector-valued basis of H_div
    # \int_\Omega div(phi_i) div(phi_j) \dx
    DivPhiDivPhi = assemble(inner(Div(y, dim), Div(q, dim)) * dx(domain=mesh))
    # \int_\Omega A^{-1} phi_i \cdot phi_j \dx
    PhiPhi = assemble(inner(invA * y, q) * dx(domain=mesh))
    # Define vectors included into the RHS
    RhsDivPhi = assemble(inner(-f_bar, Div(q, dim)) * dx(domain=mesh))
    RhsPhi = assemble(inner(Grad(v, dim), q) * dx(domain=mesh))

    #print "DivDiv = ", DivPhiDivPhi.array()
    #print "PhiPhi = ", PhiPhi.array()

    #print "RhsDivPhi = ", RhsDivPhi.array()
    #print "RhsPhi = ", RhsPhi.array()

    return DivPhiDivPhi, PhiPhi, RhsDivPhi, RhsPhi

def majorant_nd(v, y, H_div, V, VV, f, A, invA, lambda_1, a, lmbd,
                error,
                mesh, dim, C_FD, C_Ntr,
                problem_params):
    """
    :param v: approximate solution
    :param y: flux
    :param H_div: funtional space for the flux function
    :param f: right-hand side function (modified RHS)
    :param A: diffusion matrix
    :param lambda_1: minimal eigenvalue of diffusion matrix
    :param a: convection function
    :param lmbd: reaction function
    :param error: error value
    :param mesh: mesh
    :param dim: geometrical problem dimension
    :param C_FD: Freidrichs constant of the computational domain
    :param MAJORANT_OPTIMIZE: test_parameter on weather to optimize majorant of not
    :return maj: majorant value
            m_d, m_f_one_minus_mu_opt: majorant components
            beta: optimal parameter
            var_m_d, var_m_f_w_opt: variational form of the majorant (to use them later in construction of error estimators)
            maj, y, beta, m_d, m_f, var_m_d, var_m_f_w_opt, majorant_reconstruction_time, majorant_minimization_time
    """

    tic()

    # Initialize value
    beta = 1.0
    #for i in range(2):
    f_bar = f - lmbd * v - inner(a, Grad(v, dim))
    #f_bar = f
    maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = calculate_majorant_bar_mu_opt(v, y, beta, f_bar, A, invA, lambda_1, lmbd, mesh, dim, C_FD)
    i_eff = sqrt(maj / error)
    print " "
    print '%------------------------------------------------------------------------------------%'
    print "% Majorant before optimization"
    print '%------------------------------------------------------------------------------------%'
    print " "

    output_optimization_results(0, maj, m_d, m_f, i_eff, beta, error)

    majorant_reconstruction_time = toc()


    if problem_params["MAJORANT_OPTIMIZE"]:
        #----------------------------------------------------------------------------#
        # Optimization algorithm
        #----------------------------------------------------------------------------#
        print " "
        print "%-----------------------"
        print "% optimization "
        print "%-----------------------"
        print " "

        tic()
        S, K, z, g = get_matrices_of_optimization_problem_bar(H_div, v, f_bar, invA, mesh, dim)

        y = Function(H_div)
        Y = y.vector()

        OPT_ITER = problem_params["majorant_optimization_iterations"]
        # Execute iterative process to optimize majorant with respect to beta and flux
        for k in range(0, problem_params["majorant_optimization_iterations"]):
            # Solve system with respect to Y
            #yMatrix = (C_FD ** 2) / lambda_1 * S + beta * K
            #yRhs = sum((C_FD ** 2) / lambda_1 * z, beta * g)
            #solve(yMatrix, Y, yRhs)

            solve((C_FD ** 2) / lambda_1 * S + beta * K, Y, (C_FD ** 2) / lambda_1 * z + beta * g)

            y.vector = Y
            '''

            YtSY = assemble(inner(Div(y, dim), Div(y, dim)) * dx(domain=mesh))
            YtSz2 = assemble(2 * inner(Div(y, dim), f) * dx(domain=mesh))
            FF = assemble(inner(f, f) * dx(domain=mesh))

            print "Y^T*S*Y", YtSY
            print "2*Y^T*z", YtSz2
            print "FF", FF
            print "\| div y + f \|^2", YtSY + YtSz2 + FF

            YtY = assemble(inner(y, y) * dx(domain=mesh))
            YtSz2 = assemble(2 * inner(Grad(v, dim), y) * dx(domain=mesh))
            GradVGradV = assemble(inner(Grad(v, dim), Grad(v, dim)) * dx(domain=mesh))

            print "YtY", YtY
            print "YtSz2", YtSz2
            print "GradVGradV", GradVGradV
            print "\| y - grad v \|^2", YtY - YtSz2 + GradVGradV

            print "\| v \|", (norm(v, 'L2'))**2
            print "\| grad v \|", (norm(v, 'H1'))**2
            print "\| div y \|", (norm(project(div(y), V), 'L2'))**2
            print "\| f \|", (norm(project(f, V), 'L2'))**2
            '''
            # Calculate majorant
            maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = calculate_majorant_bar_mu_opt(v, y, beta, f_bar, A, invA, lambda_1, lmbd, mesh, dim, C_FD)
            i_eff = sqrt(maj / error)

            output_optimization_results(k + 1, maj, m_d, m_f, i_eff, beta, error)

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
    #dofmap = V_exact
    #mesh_coordimates = dofmap.
    # Define integration scheme order and assign the integrator
    scheme_order = 4
    gauss = integration.SpaceIntegrator(scheme_order, dim)

    # Iterate over cells of the mesh and obtain local value of the majorant
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

def error_norm(u, ue, lmbd, A, invA, a, func_a, V, Ve, mesh, dim):
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
    delta = abs(lmbd - 0.5 * Div(func_a, dim))
    var_delta_e = delta * inner(e, e)
    var_lambda_e = lmbd * inner(e, e)
    var_diva_e = (- 0.5 * Div(func_a, dim)) * inner(e, e)

    var_a_e = inner(invA * a * e, a * e)

    # Assembling variational form of the error
    grad_e = assemble(var_grad_e * dx(domain=mesh))
    delta_e = assemble(var_delta_e * dx(domain=mesh))
    lambda_e = assemble(var_lambda_e * dx(domain=mesh))
    diva_e = assemble(var_diva_e * dx(domain=mesh))
    a_e = assemble(var_a_e * dx(domain=mesh))

    # Calculate L2 norm
    l2_e = assemble(inner(e, e) * dx(domain=mesh))
    # Calculate Linf norm
    #e_func = project(e, Ve, form_compiler_parameters={'quadrature_degree': 4})
    #linf_e = norm(e_func.vector(), 'linf')

    # Calculate Linf based on the nodal values
    u_exact_v = interpolate(ue, V)
    linf_e = abs(u_exact_v.vector().array() - u.vector().array()).max()


    print '%------------------------------------------------------------------------------------%'
    print '% Error '
    print '%------------------------------------------------------------------------------------%'
    print "\| grad e \|^2_A                 = %8.2e" % grad_e
    print "\| e \|^2                        = %8.2e" % l2_e
    print "\| e \|^2_inf                    = %8.2e" % linf_e
    print "\| (lmbd - 0.5 div a)^0.5 e \|^2 = %8.2e" % delta_e
    print "\| lmbd^0.5 e \|^2               = %8.2e" % lambda_e
    print "\| a e \|^2_{A^{-1}}             = %8.2e" % a_e

    '''
    print "\| grad e \|^2_A + \| (lmbd - 0.5 div a)^0.5 e \|^2          = %8.2e" \
          % (grad_e + delta_e)
    print "\| grad e \|^2_A + \| lmbd^0.5 e \|^2 + \| a e \|^2_{A^{-1}} = %8.2e" \
          % (grad_e + lambda_e + a_e)
    '''
    return grad_e, l2_e, linf_e, delta_e, lambda_e, a_e, var_grad_e, var_delta_e, var_lambda_e, var_a_e

def minorant(v_ref, v, mesh, boundary_facets, ds, V, Vh, uD, uN, uD_boundary, f, A, lmbd, conv, dim, test_params):

    # Define variational problem
    w = TrialFunction(Vh)
    mu = TestFunction(Vh)

    def a(u, v): return inner(A * Grad(u, dim), Grad(v, dim)) + lmbd * inner(u, v)
    def b(u, v): return inner(A * Grad(u, dim), Grad(v, dim)) + lmbd * inner(u, v) + inner(inner(conv, Grad(u, dim)), v)
    def l(v):    return f * v
    def l_N(v):  return uN * v

    # Define variational forms dependent on the problem
    '''
    if test_params["pde_tag"] == "reaction-diffusion-pde" or test_params["pde_tag"] == "diffusion-pde":

        a_var = a(w, mu) * dx(domain=mesh)
        L_var = (l(mu)
                 - inner(A * Grad(v, dim), Grad(mu, dim))
                 - lmbd * inner(v, mu)) \
                * dx(domain=mesh) \
                + l_N(mu) * ds(1)
                #+ l_N(mu) * ds(domain=boundry_factes)

    elif test_params["pde_tag"] == "reaction-convection-diffusion-pde" or test_params["pde_tag"] == "convection-diffusion-pde":

        a_var = b(w, mu) * dx(domain=mesh)
        L_var = (l(mu)
                 - inner(A * Grad(v, dim), Grad(mu, dim))
                 - lmbd * inner(v, mu)
                 - inner(inner(conv, Grad(v, dim)), mu)) \
                * dx(domain=mesh) \
                + l_N(mu) * ds(1)
    '''
    if test_params["solution_tag"] == "predefined-solution":


        a_var = a(w, mu) * dx(domain=mesh)
        L_var = (l(mu) - b(v, mu)) * dx(domain=mesh) \
                + l_N(mu) * ds(1)

        # Define boundary condition and solve the minimization problem with respect to w
        w = Function(Vh)
        # Dirichlet function = 0.0 since we are looking for w =  v_ref - v
        if test_params.test_num == 45:
            bc = [DirichletBC(V, uD[0], boundary_facets, 0),  # top, right, left
                  DirichletBC(V, uD[1], boundary_facets, 1)]  # bottom
        else:
            bc = DirichletBC(Vh, Constant(0.0), uD_boundary)
        solve(a_var == L_var, w, bc)

    elif test_params["solution_tag"] == "reference-solution":
        w_ = abs(v_ref - v)
        #w_ = (v_ref - v)
        w = project(w_, Vh)

    # Define variational form for minorant: 2 * F(u, w) - a(w, w), where
    # F(u, w) = l(w) - b(v, w) // "reaction-convection-diffusion-pde", "convection-diffusion-pde"
    min_var = (2 * (l(w) - b(v, w)) - a(w, w)) * dx(domain=mesh)
    min = assemble(min_var)

    return min, w

def output_result_error_and_majorant(error, majorant, i_eff_maj):
    print ' '
    print '%------------------------------------------------------------------------------------%'
    print '% Result majorant '
    print '%------------------------------------------------------------------------------------%'
    print '||| e |||  = %8.4e ' % (error)
    print 'maj        = %8.4e' % (majorant)
    print 'i_eff(maj) = %.4f' % (i_eff_maj)
    print " "

def output_result_minorant(error, minorant, i_eff_min):
    print ' '
    print '%------------------------------------------------------------------------------------%'
    print '% Result minorant '
    print '%------------------------------------------------------------------------------------%'
    print '||| e |||   = %8.4e ' % (error)
    print 'min         = %8.4e' % (minorant)
    print 'i_eff(min)  = %.4f' % (i_eff_min)
    print " "


def bulk_marking(mesh, distr, theta):
    marking = CellFunction("bool", mesh)
    distr_0 = sorted(distr, reverse=True)[int(len(distr) * theta)]
    for cell in cells(mesh):
        marking[cell] = distr[cell.index()] > distr_0

    return marking
'''
def bulk_marking(mesh, distr, theta):
    # Sort eta_T in decreasing order and keep track of the cell numbers
    indices = distr.argsort()[::-1]
    sorted = distr[indices]
    
    # Compute sum and fraction of indicators
    total = sum(sorted)
    fraction = theta*total
    
    # Define cell function to hold markers
    markers = CellFunction("bool", mesh, False)
    
    # Iterate over the cells
    v = 0.0
    for i in indices:
        # Stop if we have marked enough
        if v >= fraction:
            break
        # Otherwise
        markers[i] = True
        v += sorted[i]
    
    return markers
'''

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

def solve_dual_problem(V_star, f, u0, boundary, u, u_e, dim):

    # Define dual variational problem
    z = TrialFunction(V_star)
    w = TestFunction(V_star)

    # Define the system
    a_star = inner(Grad(z, dim), Grad(w, dim)) * dx(domain=mesh)
    grad_e = Grad(u_e, dim) - Grad(u, dim)
    norm_grad_e = sqrt(assemble(inner(grad_e, grad_e) * dx(domain = mesh)) )
    L_star = J(w, grad_e, norm_grad_e, dim)

    z = Function(V_star)
    bc_star = DirichletBC(V_star, u0, boundary)
    solve(a_star == L_star, z, bc_star)

    z_exact = project( (u - u_e) / norm_grad_e, V_star)
    z_error = assemble((z - z_exact)**2 * dx(domain=mesh))

    print "z error", z_error

    return z

def compare_error_indicators(mesh, V, V_star, f, u0, boundary, u, u_e, grad_u_e, y, beta, test_num, tag, dim):

    z = solve_dual_problem(V_star, f, u0, boundary, u, u_e)
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

def solve_dual_problem(mesh, V_star, u0, boundary, e, dim, norm_grad_e):

    # Define dual variational problem
    z = TrialFunction(V_star)
    w = TestFunction(V_star)

    # Define the adjoint right-hand side
    a_star = inner(Grad(z, dim), Grad(w, dim)) * dx(domain=mesh)

    # Define the adjoint left-hand side
    L_star = J_energy_norm_error(w, Grad(e, dim), norm_grad_e, dim) * dx(domain=mesh)

    z = Function(V_star)
    bc_star = DirichletBC(V_star, u0, boundary)
    solve(a_star == L_star, z, bc_star)

    z_exact = project( e / norm_grad_e, V_star)

    z_error = assemble((z - z_exact)**2 * dx(domain=mesh))
    print "z error", z_error

    return z
'''
