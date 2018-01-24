from dolfin.cpp.common import tic, toc
from dolfin.cpp.function import near
from dolfin.cpp.la import list_linear_solver_methods, list_krylov_solver_preconditioners

__author__ = 'svetlana'

from dolfin import *
from dolfin.cpp.mesh import cells, BoundaryMesh, CellFunction, SubMesh, Facet, vertices
import postprocess, integration
import numpy
from problem import Grad, Div, D_t, Laplas, NablaGrad, neumann_bc_marker
import math


def output_optimization_results(iter, maj, m_d, m_f, i_eff, beta, error):
    print 'opt cycle # %d' % iter
    print 'maj^2   = %8.2e, m^2_d = %8.2e, m^2_f = %8.2e, i_eff = %.4f, beta = %.4f' \
          % (maj, m_d, m_f, i_eff, beta)
    print '[e]^2 = %8.2e' % (error)

#def update_majorant_components(u, y, f, c_H, mesh, dim):
def update_majorant_components(r_d, r_f, invA, eps, mesh):

    m_d = assemble(inner(invA * r_d, r_d) * dx(domain=mesh))
    m_f = assemble(inner(r_f, r_f) * dx(domain=mesh))

    return m_d, m_f

def calculate_majorant_bar_mu_opt(u, y, beta, f_bar, A, invA, min_eig_A, lmbd, mesh, dim, C_FD):
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
    mu_opt = (C_FD ** 2) * (1.0 + beta) * lmbd / (beta * min_eig_A + (C_FD ** 2) * (1.0 + beta) * lmbd)
    #mu_opt = 1
    w_opt = (C_FD ** 2) * (1.0 + beta) / (beta * min_eig_A + (C_FD ** 2) * (1.0 + beta) * lmbd) # this is the function in front of inner(r_f, r_f)

    # Define residuals
    r_d = y - A * NablaGrad(u, dim)
    r_f = (Div(y, dim) + f_bar)

    #r_d = y - eps * A * Grad(u, dim)
    #r_f = Div(y, dim) + f_bar

    # Define variational forms
    var_m_d = inner(invA * r_d, r_d)
    var_m_f_w_opt = w_opt * inner(r_f, r_f)
    var_m_f_one_minus_mu_opt = ((1 - mu_opt) ** 2) * inner(r_f, r_f) # only for the calculation of optimal beta

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
    beta = C_FD * sqrt(m_f_one_minus_mu_opt / m_d / min_eig_A)

    return maj, m_d, m_f_one_minus_mu_opt, beta, var_m_d, var_m_f_w_opt

def update_majorant_II_components(u, w, y, f, u0, mesh, dim, v_deg, W, Ve, T):

    r_d = (y - Grad(u, dim) + Grad(w, dim))
    r_f = (Div(y, dim) + f - D_t(u, dim) - D_t(w, dim))
    L = D_t(u, dim) * w + inner(Grad(u, dim), Grad(w, dim)) - f * w
    w_T, T_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, w, T, dim, v_deg+1)
    w_0, O_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, w, 0, dim, v_deg+1)

    u0_0, O_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, interpolate(u0, Ve), 0, dim, v_deg) # varphi
    u_0, O_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, interpolate(u, Ve), 0, dim, v_deg)
    r_0 = u0_0 - u_0
    l = inner(r_0, r_0) - 2 * w_0 * r_0
    w_T_form = inner(w_T, w_T)


    m_d = assemble(inner(r_d, r_d) * dx(domain=mesh))
    m_f = assemble(inner(r_f, r_f) * dx(domain=mesh))
    m_L = assemble((D_t(u, dim) * w + inner(Grad(u, dim), Grad(w, dim)) - f * w) * dx(domain=mesh))
    m_T = assemble(w_T_form * dx(domain=T_mesh))
    m_l = assemble(l * dx(domain=O_mesh))

    return m_d, m_f, m_l, m_L, m_T

def calculate_majorant(m_d, m_f, m_b, eps, C_FD):

    if m_f <= DOLFIN_EPS:
        beta = 1.0
        maj = m_b + m_d
    else:
        if m_d <= DOLFIN_EPS:
            beta = 1.0
            maj = m_b + (C_FD ** 2) / eps * m_f
        else:
            beta = C_FD * sqrt(m_f / m_d / eps)
            maj = m_b + (1.0 + beta) * m_d + (1.0 + 1.0 / beta) * (C_FD ** 2) / eps * m_f

    return maj, beta

def calculate_majorant_II(m_d, m_f, m_b, C_FD):

    if m_f <= DOLFIN_EPS:
        beta = 1.0
        maj = m_b + m_d
    else:
       if m_d <= DOLFIN_EPS:
           beta = 1.0
           maj = m_b + (C_FD ** 2) * m_f
       else:
           beta = C_FD * sqrt(m_f / m_d)
           maj = m_b + (1.0 + beta) * m_d + (1.0 + 1.0 / beta) * (C_FD ** 2) * m_f

    return maj, beta

def calculate_majorant_II(m_d, m_f, l, L, m_T, C_FD, gamma):

    if m_f <= DOLFIN_EPS:
        beta = 1.0
        maj_II = gamma * m_T + l + 2 * L + m_d
    else:
       if m_d <= DOLFIN_EPS:
           beta = 1.0
           maj_II = gamma * m_T + l + 2 * L + (C_FD ** 2) * m_f
       else:
           beta = C_FD * sqrt(m_f / m_d)
           maj_II = gamma * m_T + l + 2 * L + (1.0 + beta) * m_d + (1.0 + 1.0 / beta) * (C_FD ** 2) * m_f

    return maj_II, beta


#def get_matrices_of_optimization_problem(H_div, u, f, c_H, mesh, dim):
def get_matrices_of_optimization_problem(H_div, u, f_bar, invA, mesh, dim):

    # Define variational problem
    y = TrialFunction(H_div)
    q = TestFunction(H_div)

    # Define system of linear equation to find the majorant
    S = assemble(inner(Div(y, dim), Div(q, dim)) * dx(domain=mesh))
    K = assemble(inner(invA * y, q) * dx(domain=mesh))
    N = assemble(inner(y, q) * ds(neumann_bc_marker))

    #z = assemble(inner(- (f - c_H * D_t(u, dim)), Div(q, dim)) * dx(domain=mesh))
    z = assemble(inner(-f_bar, Div(q, dim)) * dx(domain=mesh))
    g = assemble(inner(NablaGrad(u, dim), q) * dx(domain=mesh))

    return S, K, z, g

def get_matrices_of_optimization_problem_II(W, u, y, f, u0, T, mesh, dim, v_deg):


    # Define variational problem
    w = TrialFunction(W)
    mu = TestFunction(W)

    #w_T, T_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, w, T, dim, v_deg + 1)
    #mu_T, T_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, mu, T, dim, v_deg + 1)

    #mu_0, O_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, mu, 0, dim, v_deg + 1)
    #u0_0, O_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, interpolate(u0, W), 0, dim, v_deg + 1)


    # Define system of linear equation to find the majorant
    S = assemble(inner(Grad(w, dim), Grad(mu, dim)) * dx(domain=mesh))
    K = assemble(inner(D_t(w, dim), D_t(mu, dim)) * dx(domain=mesh))
    #F = assemble(inner(w_T, mu_T)) * dx(domain=T_mesh)

    L = assemble((D_t(u, dim) * mu + inner(Grad(u, dim), Grad(mu, dim)) - f * mu) * dx(domain=mesh))
    z = assemble((inner(y - Grad(u, dim), Grad(mu, dim))) * dx(domain=mesh))
    g = assemble((f - D_t(u, dim) + Div(y, dim)) * D_t(mu, dim) * dx(domain=mesh))
    #I = assemble(u0_0 * mu_0 * dx(domain=O_mesh))

    return S, K, L, z, g

def calculate_CF_of_domain(domain, dim):

    if domain == "l-shape-domain" and dim == 3:
        height = 1
        width = 2
        length = 2
        C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height**2 + 1.0 / width**2 + 1.0 / length**2)

    elif domain == "1times1-minus-0times0" and dim == 3:
        height = 2
        width = 2
        length = 2
        C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height**2 + 1.0 / width**2 + 1.0 / length**2)

    elif domain == "unit-domain":

        height = 1
        width = 1
        length = 1

        if dim == 3:
            C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height**2 + 1.0 / width**2)
        elif dim == 2:
            C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / length**2)

    elif domain == "l-shape-domain" and dim == 2:
        width = 2
        length = 2
        C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / width**2 + 1.0 / length**2)


    print "C_FD = 1 / pi / (h_1^(-2) + ... + h_n^(-2)) = ", C_FD

    return C_FD


def majorant_nd(u, Ve, y, H_div, f, A, invA, min_eig_A, c_H, eps, lmbd, a, u0, error, mesh, C_FD, dim, test_params):

    tic()

    u0_0, O_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, interpolate(u0, Ve), 0, dim, test_params["v_approx_order"],)
    u_0, O_mesh = get_2d_slice_of_3d_function_on_Oz(mesh, interpolate(u, Ve), 0, dim, test_params["v_approx_order"],)
    r_b = abs(u0_0 - u_0)
    m_b = assemble(inner(r_b, r_b) * dx(domain=O_mesh))

    # Define residuals
    beta = 1.0
    f_bar = f - c_H * D_t(u, dim) - lmbd * u - inner(a, NablaGrad(u, dim))
    maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = calculate_majorant_bar_mu_opt(u, y, beta, f_bar,
                                                                                A, invA, min_eig_A,
                                                                                lmbd, mesh, dim, C_FD)

    i_eff = sqrt(maj / error)

    print " "
    print '%------------------------------------------------------------------------------------%'
    print "% Majorant before optimization"
    print '%------------------------------------------------------------------------------------%'
    print " "

    majorant_reconstruction_time = toc()
    output_optimization_results(-1, maj, m_d, m_f, i_eff, beta, error)

    if test_params['MAJORANT_OPTIMIZE']:
        #----------------------------------------------------------------------------#
        # Optimization algorithm
        #----------------------------------------------------------------------------#
        print " "
        print "%-----------------------"
        print "% optimization "
        print "%-----------------------"
        print " "

        tic()
        #S, K, z, g = get_matrices_of_optimization_problem(H_div, u, f, c_H, mesh, dim)
        S, K, z, g = get_matrices_of_optimization_problem(H_div, u, f_bar, invA, mesh, dim)

        y = Function(H_div)
        Y = y.vector()

        # Execute iterative process to optimize majorant with respect to beta and flux
        for k in range(1, test_params['majorant_optimization_iterations']):

            #list_linear_solver_methods()
            #list_krylov_solver_preconditioners()

            # Solve system with respect to Y
            solve((C_FD ** 2) / min_eig_A * S + beta * K, Y, (C_FD ** 2) / min_eig_A * z + beta * g) # lu 3.4 - 3.5
            #solve((C_FD ** 2) * S + beta * K, Y, (C_FD ** 2) * z + beta * g, "gmres", "jacobi") # 3.4
            #solve((C_FD ** 2) * S + beta * K, Y, (C_FD ** 2) * z + beta * g, "cg", "ilu") # 3.4
            #solve(C_FD * S + beta * K, Y, C_FD * z + beta * g, "cg", "jacobi")
            #solve(C_FD * S + beta * K, Y, C_FD * z + beta * g, "gmres", "hypre_amg") # 3.4 - 3.5

            y.vector = Y
            # Update majorant
            maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = calculate_majorant_bar_mu_opt(u, y, beta, f_bar, A, invA,
                                                                                        min_eig_A,
                                                                                        lmbd, mesh, dim, C_FD)
            i_eff = sqrt(maj / error)

            output_optimization_results(k, maj, m_d, m_f, i_eff, beta, error)

        majorant_minimization_time = toc()
    else:

        majorant_minimization_time = 0.0

    return maj, y, beta, m_d, m_f, var_m_d, var_m_f_w_opt, majorant_reconstruction_time, majorant_minimization_time



def majorant_II_nd(u,  Ve, w, W, y, f, u0, error_II, gamma, T, mesh, C_FD, dim, v_deg, MAJORANT_OPTIMIZE):

    tic()

    m_d, m_f, l, L, m_T = update_majorant_II_components(u, w, y, f, u0, mesh, dim, v_deg, W, Ve, T)
    maj_II, beta = calculate_majorant_II(m_d, m_f, l, L, m_T, C_FD, gamma)
    i_eff = sqrt(maj_II / error_II)

    majorant_reconstruction_time = toc()
    output_optimization_results(-1, maj_II, m_d, m_f, i_eff, beta, error_II)

    if MAJORANT_OPTIMIZE:
        #----------------------------------------------------------------------------#
        # Optimization algorithm
        #----------------------------------------------------------------------------#
        tic()
        S, K, L, z, g = get_matrices_of_optimization_problem_II(W, u, y, f, u0, T, mesh, dim, v_deg)

        w = Function(W)
        W = w.vector()

        OPT_ITER = 4
        # Execute iterative process to optimize majorant with respect to beta and flux
        for k in range(1, OPT_ITER):

            #list_linear_solver_methods()
            #list_krylov_solver_preconditioners()

            # Solve system with respect to Y
            solve((1 + beta) * S - (C_FD ** 2) * beta / (1 + beta) * K, W,
                  - L - (1 + beta) * z - (C_FD ** 2) * beta / (1 + beta) * g) # lu 3.4 - 3.5
            # Calculate majorant
            m_d, m_f, l, L, m_T = update_majorant_II_components(u, w, y, f, u0, mesh, dim, v_deg, W, Ve, T)
            maj_II, beta = calculate_majorant_II(m_d, m_f, l, L, m_T, C_FD, gamma)
            i_eff = sqrt(maj_II / error_II)

            output_optimization_results(k, maj_II, m_d, m_f, i_eff, beta, error_II)

        majorant_minimization_time = toc()
    else:

        majorant_minimization_time = 0.0


    return maj_II, beta, majorant_reconstruction_time, majorant_minimization_time

def error_majorant_distribution_nd(mesh, dim, V_exact,
                                   var_grad_e, var_delta_e, var_m_d, var_m_f_w_opt,
                                   beta, C_FD):
    cell_num = mesh.num_cells()
    ed_distr = postprocess.allocate_array(cell_num)
    delta_e_distr = postprocess.allocate_array(cell_num)
    md_distr = postprocess.allocate_array(cell_num)
    maj_distr = postprocess.allocate_array(cell_num)

    # Project UFL forms on the high-order functional space to obtain functions
    ed = project((var_grad_e), V_exact)
    delta_e = project(var_delta_e, V_exact)
    md = project((var_m_d), V_exact)
    mf_wopt = project(var_m_f_w_opt, V_exact)

    scheme_order = 4
    gauss = integration.SpaceIntegrator(scheme_order, dim)

    for c in cells(mesh):
        # Obtaining the coordinates of the vertices in the cell
        verts = c.get_vertex_coordinates().reshape(dim + 1, dim)
        x_n_transpose = verts.T

        # Calculating the area of the cell
        matrix = postprocess.allocate_array_2d(dim + 1, dim + 1)
        matrix[0:dim, 0:dim + 1] = x_n_transpose
        matrix[dim, 0:dim + 1] = numpy.array([1.0 for i in range(dim + 1)])
        # print "matrix = ", matrix
        meas = abs(numpy.linalg.det(matrix))

        # Integrating over the cell
        ed_distr[c.index()] = gauss.integrate(ed, meas, x_n_transpose)
        delta_e_distr[c.index()] = gauss.integrate(delta_e, meas, x_n_transpose)
        md_distr[c.index()] = gauss.integrate(md, meas, x_n_transpose)
        #maj_distr[c.index()] = gauss.integrate(mdf, meas, x_n_transpose)

    #print 'num md_cells', md_distr
    #print 'ed_cells', ed_distr
    #print 'delta_e_cells', delta_e_distr

    #print 'maj_cells', maj_distr
    #print 'e_cells', e_distr

    print '\nmd = Sum_K md_K = %8.2e' % sum(md_distr)
    print 'ed = Sum_K ed_K = %8.2e' % sum(ed_distr)

    return ed_distr, md_distr, maj_distr

def get_2d_slice_of_3d_function_on_Oz(mesh, u, T, dim, v_degree):

    # create a boundary mesh
    bmesh = BoundaryMesh(mesh, "exterior")
    cell_func = CellFunction('size_t', bmesh, 0)
    coordinates = bmesh.coordinates()
    z_indeces = postprocess.allocate_array(dim)

    for cell in cells(bmesh):
        indx = 0
        for vertex in vertices(cell):
            z_indeces[indx] = coordinates[vertex.index()][dim-1]
            indx += 1
            #print "Vertex index with coordinates", vertex.index(), coordinates[vertex.index()][dim-1]
        if (dim == 3 and near(z_indeces[0], T) and near(z_indeces[1], T) and near(z_indeces[2], T)) or (dim == 2 and near(z_indeces[0], T) and near(z_indeces[1], T)) :
            #print "right cell", cell.index()
            cell_func[cell] = 1

    submesh = SubMesh(bmesh, cell_func, 1)

    # create a FunctionSpace on the submesh-
    Vs = FunctionSpace(submesh, "Lagrange", v_degree)
    us = interpolate(u, Vs)

    return us, submesh

def error_norm(u, ue, A, lmbd, a, f, c_H, v_degree, mesh, T, dim, V, Ve):

    """
    :param u: approximate solution
    :param ue: exact solution
    :param A:
    :param lmbd:
    :param lambd:

    :param Ve: functional space of exact solution
    :return: L2 error-norm between u and ue
    """
    u_ve = interpolate(u, Ve)
    u_exact_ve = interpolate(ue, Ve)
    e = abs(u_ve - u_exact_ve)
    #e = (u_ve - u_exact_ve)

    res = Div(A * NablaGrad(u, dim), dim) \
          + f \
          - c_H * D_t(u, dim) \
          - lmbd * u \
          - inner(a, NablaGrad(u, dim))

    var_e        = inner(e, e)
    var_delta_e  = inner((lmbd - 0.5 * Div(a, dim)) * e, e)
    var_lambda_e = lmbd * inner(e, e)
    var_grad_e   = inner(A * NablaGrad(e, dim), NablaGrad(e, dim))
    var_e_t      = inner(c_H * D_t(e, dim), c_H * D_t(e, dim))
    var_laplas_e = inner(Div(A * NablaGrad(e, dim), dim), Div(A * NablaGrad(e, dim), dim))
    var_e_id     = inner(res, res)

    val_e        = assemble(var_e * dx(domain=mesh))
    val_delta_e  = assemble(var_delta_e * dx(domain=mesh))
    val_lmbd_e   = assemble(var_lambda_e * dx(domain=mesh))
    val_grad_e   = assemble(var_grad_e * dx(domain=mesh))
    val_e_t      = assemble(var_e_t * dx(domain=mesh))
    val_laplas_e = assemble(var_laplas_e * dx(domain=mesh))

    val_e_id     = assemble(var_e_id * dx(domain=mesh))

    u_T, mesh_T  = get_2d_slice_of_3d_function_on_Oz(mesh, u_ve, T, dim, v_degree)
    ue_T, mesh_T = get_2d_slice_of_3d_function_on_Oz(mesh, u_exact_ve, T, dim, v_degree)
    var_e_T      = abs(ue_T - u_T)
    val_e_T      = assemble(inner(var_e_T, var_e_T) * dx(domain=mesh_T))

    u_0, mesh_0  = get_2d_slice_of_3d_function_on_Oz(mesh, u_ve, 0.0, dim, v_degree)
    ue_0, mesh_0 = get_2d_slice_of_3d_function_on_Oz(mesh, u_exact_ve, 0.0, dim, v_degree)
    var_e_0 = abs(ue_0 - u_0)
    val_e_0 = assemble(inner(A * var_e_0, var_e_0) * dx(domain=mesh_0))

    print '%------------------------------------------------------------------------------------%'
    print '% Error '
    print '%------------------------------------------------------------------------------------%\n'
    print "\| grad_x e \|^2_A               = %8.2e" % val_grad_e
    print "\| e \|^2                        = %8.2e" % val_e
    print "\| (lmbd - 0.5 div a)^0.5 e \|^2 = %8.2e" % val_delta_e
    print "\| lmbd^0.5 e \|^2               = %8.2e\n" % val_lmbd_e

    print "\| laplas e \|^2                 = %8.2e" % val_laplas_e
    print "\| e_t \|^2                      = %8.2e" % val_e_t
    print "\| r_v \|^2                      = %8.2e" % val_e_id

    print "\| e \|^2_T                      = %8.2e" % val_e_T
    print "\| e \|^2_0                      = %8.2e\n" % val_e_0
    #print "\| grad_x e \|^2_T               = %8.2e" % val_grad_e_T
    #print "\| grad_x e \|^2_0               = %8.2e\n" % val_grad_e_0

    print '%------------------------------------------------------------------------------------%'
    print '% Error identity '
    print '%------------------------------------------------------------------------------------%\n'

    print "id = \| e \|^2_T + \| r_v \|^2               = %8.2e" % (val_e_T + val_e_id)
    print "\| laplas e \|^2 + \| e_t \|^2 + \| e \|^2_0 = %8.2e\n" % (val_e_0 + val_laplas_e + val_e_t + val_delta_e)

    return e, val_e, val_grad_e, val_delta_e, val_e_T, val_laplas_e + val_e_t, val_e_id, \
           var_e, var_delta_e, var_grad_e

def minorant(u, mesh, Vh, u0, u0_boundary, f, dim, error):

    # Define variational problem
    w = TrialFunction(Vh)
    mu = TestFunction(Vh)

    a = inner(NablaGrad(w, dim), NablaGrad(mu, dim)) * dx(domain=mesh)
    L = (f * mu - inner(NablaGrad(u, dim), Grad(mu, dim))) * dx(domain=mesh)

    w = Function(Vh)

    # Define boundary condition
    bc = DirichletBC(Vh, u0, u0_boundary)
    solve(a == L, w, bc)

    min_var = (2 * ( f * w - inner(Grad(u, dim), Grad(w, dim))) - inner(Grad(w, dim), Grad(w, dim))) * dx(domain=mesh)
    min = assemble(min_var)

    print "error = %.4e , min = %.4e, i_eff_min = %.4e" % (error, min, min/error)

    return min, w

def output_result_error_and_majorant(error, majorant, i_eff_maj):

    print '\n%------------------------------------------------------------------------------------%'
    print '% Majorant '
    print '%------------------------------------------------------------------------------------%\n'

    print '[e]       = %8.4e ' % (error)
    print 'maj       = %8.4e' % (majorant)
    print 'i_eff_maj = %.4f\n' % (i_eff_maj)

def output_result_minorant(minorant, i_eff_min):
    print 'minorant  = %8.4e' % (minorant)
    print 'i_eff_min = %.4f' % (i_eff_min)

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


def get_indicators_CG0(e_form, norm_grad_e, V, V_star, f, u0, u0_boundary, u, u_e, mesh, dim):

    e_form = u_e - u
    # Contruct the solution of the adjoint problem
    z = solve_dual_problem(mesh, V_star, u0, u0_boundary, e_form, dim, norm_grad_e)
    Pz = project(z, V)

    # Get parameters of the mesh
    h = mesh.hmax()
    n = FacetNormal(mesh)

    # Construct the set of
    DG0 = FunctionSpace(mesh, "DG", 0)
    w = TestFunction(DG0)

    r_h = jump(Grad(u, dim), n)
    R_h = Div(Grad(u, dim), dim) + f

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

def majorant_distribution_using_DG0(e_form, r_d, r_f, A, invA, min_eig_A, beta, C_FD, mesh, dim):

    # Define the functional space used for distribution over cells
    DG0 = FunctionSpace(mesh, "DG", 0)
    w = TestFunction(DG0)

    m_d_DG0 = Function(DG0)
    m_df_DG0 = Function(DG0)
    e_DG0 = Function(DG0)

    # Define variational forms of dual component of majorant and the whole functional
    m_d_var = w * (inner(invA * r_d, r_d)) * dx(domain=mesh)
    m_df_var = w * ((1.0 + beta) * inner(invA * r_d, r_d) + C_FD**2 / min_eig_A * (1 + 1 / beta) * inner(r_f, r_f)) * dx(domain=mesh)
    e_var = w * (inner(A * Grad(e_form, dim), Grad(e_form, dim))) * dx(domain=mesh)

    # Assemble the variation form and dumping obtained vector into function from DG0
    assemble(m_d_var, tensor=m_d_DG0.vector())
    assemble(m_df_var, tensor=m_df_DG0.vector())
    assemble(e_var, tensor=e_DG0.vector())

    # Assign distributions
    m_d_distr = m_d_DG0.vector().array()
    m_df_distr = m_df_DG0.vector().array()
    e_distr = e_DG0.vector().array()

    # Check the correctness of the localized distributions
    print "sum(maj_distr)", numpy.sum(m_df_distr)
    print "sum(m_d_distr)", numpy.sum(m_d_distr)
    print "sum(e_d_distr)", numpy.sum(e_distr)

    return e_distr, m_d_distr, m_df_distr

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
    fraction = theta * total

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

def averaged_marking(mesh, distr):
    cells_num = mesh.num_cells()
    marking = CellFunction("bool", mesh)
    distr_aver = sum(distr) / cells_num
    for c in cells(mesh):
        marking[c] = distr[c.index()] >= distr_aver
        #print distr[c.index()] >= distr_aver
    return marking
