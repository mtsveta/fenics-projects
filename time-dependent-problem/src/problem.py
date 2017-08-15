__author__ = 'svetlana'

__author__ = 'svetlana'

from dolfin import *
from numpy import info

'''
def construct_stiffness_and_mass_matrices(u, v, lmbd):

    # Components of stiffness (K) & mass (M) matrices
    a_K = inner(nabla_grad(u), nabla_grad(v)) * dx
    a_M = (1 + lmbd) * u * v * dx

    K = assemble(a_K)
    M = assemble(a_M)

    return K, M
'''

def construct_stiffness_and_mass_matrices(u, v, dim):
    # Components of stiffness (K) & mass (M) matrices
    a_K = inner(Grad(u, dim), Grad(v, dim)) * dx
    a_M = u * v * dx

    K = assemble(a_K)
    M = assemble(a_M)
    return K, M


def functional_spaces(mesh, problem_params, dim):

    # Order of the space of the exact solutions
    v_degree = problem_params["v_approx_order"]
    y_degree = problem_params["flux_approx_order"]
    v_exact_degree = problem_params["v_exact_approx_order"]

    # Define functional spaces
    V = FunctionSpace(mesh, 'Lagrange', v_degree)
    V_exact = FunctionSpace(mesh, 'Lagrange', v_exact_degree)

    if dim == 1:
        H_div = FunctionSpace(mesh, 'Lagrange', y_degree)
        VV = FunctionSpace(mesh, 'Lagrange', v_degree)
        VV_exact = FunctionSpace(mesh, 'Lagrange', v_exact_degree)

    else:
        H_div = FunctionSpace(mesh, 'RT', y_degree)
        #H_div = VectorFunctionSpace(mesh, 'Lagrange', y_degree)
        VV = VectorFunctionSpace(mesh, 'Lagrange', v_degree)
        VV_exact = VectorFunctionSpace(mesh, 'Lagrange', v_exact_degree)

    return V, VV, V_exact, VV_exact, H_div


# Output description of the problem
def output_problem_characteristics(test_num, u_expr, grad_u_expr, f_expr, A_expr, lambda_expr, a_expr, uD_expr,
                                   domain, T, dim, num_cells, num_vertices, nt, space_approx_deg, flux_approx_deg, delta):
    print " "
    print "%-------------------------"
    print "% Problem characteristics"
    print "%-------------------------"

    print "test: ", test_num
    print "domain: ", domain
    print "T: ", T
    print "f = ", f_expr
    print "lambda = ", lambda_expr
    print "A = ", A_expr
    print "a = ", a_expr
    print "uD = ", uD_expr
    print "u = ", u_expr
    print "grad_u = ", grad_u_expr

    print 'mesh parameters: nt = %d, num_cells = %d, num_vertices = %d' % (nt, num_cells, num_vertices)
    print "space func approx_deg = ", space_approx_deg
    print "flux approx_deg = ", flux_approx_deg
    print "dim", dim

def calculate_CF_of_domain(domain, dim):

    # We can initilize the length of the side of the domain big enough so that the 1/length is of a machine' zero order
    height, width, length = (1e16, 1e16, 1e16)

    if   domain == "l-shape-domain":
        if   dim == 3: height, width, length = (1.0, 2.0, 2.0)
        elif dim == 2: (width, length) = (2.0, 2.0)
    elif domain == "1times1-minus-0times0" and dim == 3: height, width, length = (2.0, 2.0, 2.0)
    elif domain == "unit-domain":
        if   dim == 3: height, width, length = (1.0, 1.0, 1.0)
        elif dim == 2: (width, length) = (1.0, 1.0)
        elif dim == 1: length = 1.0
    elif domain == "pi-shape-domain" and dim == 2: width, length = (2.0, 2.0)
    elif domain == "circle-domain": width, length = (2.0, 2.0)  # automatically 2d

    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height ** 2 + 1.0 / width ** 2 + 1.0 / length ** 2)
    print "C_FD = 1 / pi / sqrt(h_1^(-2) + ... + h_n^(-2)) = ", C_FD

    return C_FD

def Grad(f, dim):
    """
    :param f: function
    :param dim: dimension of the problem
    :return: gradient of the function depending on dim
    """
    if dim > 1:
        div_f = grad(f)
    else:
        # Define grad operator for the 1D case
        div_f = Dx(f, 0)
    return div_f

def Div(f, dim):
    """
    :param f: function
    :param dim: dimension of the problem
    :return: divergence of the function depending on dim
    """
    if dim > 1:
        div_f = div(f)
    else:
        # Define divergent operator for the 1D case
        div_f = Dx(f, 0)
    return div_f

def get_next_step_solution(v_k, v_k1, K, tau, M, f_k, f_k1, bc, time_discretization_tag):
    if time_discretization_tag == 'implicit':

        # IMPLICIT time-discretization scheme
        # Compose the matrix of system of linear equations
        A = M + tau * K
        # Update the right-hand side of the system of linear equation
        b = M * v_k.vector() + tau * M * f_k1.vector()

    elif time_discretization_tag == 'explicit':
        # EXPLICIT time-discretization scheme
        # Update right hand side f_{k + 1}
        # Compose the matrix of system of linear equations
        A = M
        # Update the right-hand side of the system of linear equation
        b = (M + tau * K)* v_k.vector() + tau * M * f_k.vector()


    # Solve system of linear equations
    bc.apply(A, b)
    solve(A, v_k1.vector(), b)

    return v_k1
