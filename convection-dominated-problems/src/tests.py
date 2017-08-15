__author__ = "Svetlana Matculevich <svetlana.v.matculevich@jyu.fi>"
__date__ = ""
__copyright__ = ""
__license__ = ""

import numpy
from numpy.linalg import eig
# checked
def quadratic_polynomial_solution_1d_sigma_1000():

    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    eps = 1
    A_expr = ["eps"]
    min_eig_A = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 / (1000.0 * sqrt(2 * pi)) * exp( - pow(x[0] - 0.5, 2) / (2.0 * pow(1000.0, 2)))"

    # div(A * grad u)
    divA_gradu_expr = "(" + A_expr[0] + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

# checked
def linear_trigonometric_solution_1d_rho_1minus3():

    # Solution and it's derivatives
    u_expr = "x[0]*sin(3*pi*x[0])"
    du_dx_expr = "sin(3*pi*x[0]) + x[0]*3*pi*cos(3*pi*x[0])"
    d2u_dx2_expr = "6*pi*cos(3*pi*x[0]) - 9*pow(pi, 2)*x[0]*sin(3*pi*x[0])"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    eps = 1
    A_expr = ["eps"]
    min_eig_A = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "0.001 * (x[0] + 0.001)"

    # div(A * grad u)
    divA_gradu_expr = "(" + A_expr[0] + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


# checked
def linear_trigonometric_solution_1d_rho_1():

    # Solution and it's derivatives
    u_expr = "x[0]*sin(3*pi*x[0])"
    du_dx_expr = "sin(3*pi*x[0]) + x[0]*3*pi*cos(3*pi*x[0])"
    d2u_dx2_expr = "6*pi*cos(3*pi*x[0]) - 9*pow(pi, 2)*x[0]*sin(3*pi*x[0])"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    eps = 1.0
    A_expr = ["eps"]
    min_eig_A = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 * (x[0] + 0.001)"

    # div(A * grad u)
    divA_gradu_expr = "(" + A_expr[0] + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


# checked
def linear_trigonometric_solution_1d_rho_1000():

    # Solution and it's derivatives
    u_expr = "x[0]*sin(3*pi*x[0])"
    du_dx_expr = "sin(3*pi*x[0]) + x[0]*3*pi*cos(3*pi*x[0])"
    d2u_dx2_expr = "6*pi*cos(3*pi*x[0]) - 9*pow(pi, 2)*x[0]*sin(3*pi*x[0])"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    eps = 1.0
    A_expr = ["eps"]
    min_eig_A = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1000.0 * (x[0] + 0.001)"

    # div(A * grad u)
    divA_gradu_expr = "(" + A_expr[0] + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

# checked
def quadratic_polynomial_solution_1d_sigma_1_epsilon_1eminus2():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    eps = 1e-2
    A_expr = ["eps"]
    min_eig_A = eps
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 / (1.0 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(1.0, 2)))"

    # div(A * grad u)
    divA_gradu_expr = "(" + A_expr[0] + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


# checked
def quadratic_polynomial_solution_1d_sigma_1minus2():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    eps = 1.0
    A_expr = ["eps"]
    min_eig_A = eps
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 / (0.01 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(0.01, 2)))"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"
    
    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

def quadratic_polynomial_solution_1d_lmdb_zero_a_cubic():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    eps = 1.0
    A_expr = ["eps"]
    min_eig_A = eps
    a_expr = ["-100*x[0]*x[0]*x[0]"]
    lambda_expr = "0.0"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + A_expr[0] + ")"

    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"
    
    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

# checked
def classic_example_1d_stynes_example_1_2():

    eps = 1e-5
    u_expr = "x[0] - ( exp((x[0] - 1)/(eps)) - exp(-1/(eps)) )/(1 - exp(-1/(eps)))"
    du_dx_expr = "1 - 1 / (eps) * ( exp( (x[0] - 1) / (eps) ) )/(1 - exp( - 1 / (eps) ))"
    d2u_dx2_expr = "- 1 / (eps*eps) * ( exp( (x[0] - 1) / (eps) ) )/(1 - exp( - 1 / (eps) ))"

    grad_u_expr = [du_dx_expr]
   
    A_expr = [["eps"]]
    A_inv_expr = [["1/eps"]]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["1.0"]

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, A_inv_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

# checked
def classic_example_1d_b_1_eps_1eminis2():

    eps = 0.001
    b = "1.0"
    u_expr = "(1 - exp((" + b + ") * x[0] / (eps)))/(1 - exp(" + "(" + b + ") / (eps)))"
    du_dx_expr = "-(" + b + ") / (eps) * exp((" + b + ") * x[0] / (eps))/( 1 - exp(" + "(" + b + ") / (eps)))"
    grad_u_expr = [du_dx_expr]
   
    A_expr = ["eps"]
    A_inv_expr = ["1/eps"]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = [b]

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "0.0"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, A_inv_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


# checked
def classic_example_1d_b_3_eps_1eminis2():
    # classic_example_1d_b_1_eps_1eminis2

    eps = 1e-2
    b = "3.0"
    u_expr = "(1 - exp((" + b + ") * x[0] / (eps)))/(1 - exp(" + "(" + b + ") / (eps)))"
    du_dx_expr = "-(" + b + ") / (eps) * exp((" + b + ") * x[0] / (eps))/( 1 - exp(" + "(" + b + ") / (eps)))"
    grad_u_expr = [du_dx_expr]

    A_expr = ["eps"]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = [b]

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "0.0"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def f_1_uD_polynomial_1d_lmdb_zero_a_const_eps_1minus3():

    # Input data
    f_expr = "1.0"
    eps = 1e-3
    A_expr = ["eps"]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["1.0"]

    # Dirichlet BC
    uD_expr = "0.0"
    uN_expr = "0.0"

    # Domain
    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "convectio-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def f_1_uD_zero_1d_lmdb_zero_a_cubic_eps_1minus3():

    # Input data
    eps = 1e-3
    f_expr = "1.0"
    A_expr = ["eps"]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["-(x[0] + 1.0)*(x[0] + 1.0)*(x[0] + 1.0)"]

    # Dirichlet BC
    uD_expr = "0.0"
    uN_expr = "0.0"

    # Domain
    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def f_1_uD_zero_2d_lmdb_zero_a_const_eps_1minus2():

    # Input data
    f_expr = "1.0"
    eps = 1e-2
    A_expr = [["eps",  "0"], ["0", "eps"]]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["1", "0"]

    # Dirichlet BC
    uD_expr = "0.0"
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def f_1_uD_zero_2d_lmdb_zero_a_const_eps_1minus3():

    f_expr = "1.0"
    eps = 1e-3
    A_expr = [["eps",  "0"], ["0", "eps"]]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["1", "1"]

    # Dirichlet BC
    uD_expr = "0.0"
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

def example_1_book_stynes_2d():

    f_expr = "1.0"
    eps = 1e-3
    A_expr = [["eps",  "0"], ["0", "eps"]]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["0", "0"]

    # Dirichlet BC
    uD_expr = "0"
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def example_2_eps_1eminus3book_stynes_2d():

    f_expr = "x[0]"
    eps = 1e-3
    A_expr = [["eps",  "0"], ["0", "eps"]]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["-x[0]", "0"]

    # Dirichlet BC
    uD_expr = "0"
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def example_2_eps_1eminus2book_stynes_2d():

    f_expr = "x[0]"
    eps = 1e-2
    A_expr = [["eps",  "0"], ["0", "eps"]]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["-x[0]", "0"]

    # Dirichlet BC
    uD_expr = "0"
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def example_3_book_stynes_2d():

    f_expr = "1.0"
    eps = 1e-3
    A_expr = [["eps",  "0"], ["0", "eps"]]
    min_eig_A = eps
    lambda_expr = "10.0"
    a_expr = ["x[0]*x[0] - 1", "x[1]*x[1] - 1"]

    # Dirichlet BC
    uD_expr = "0"
    uN_expr = "0.0"

    domain = "circle-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

def example_6_kleiss_tomar_2d():

    f_expr = "(x[0] - 0.1)**2 + (x[1] - 0.5)**2 <= 0.3 ? 1.0 : 0.0"
    eps = 1e-1
    A_expr = [["eps",  "0"], ["0", "eps"]]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["100", "100"]

    # Dirichlet BC
    uD_expr = ["0", "1"]
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

# checked
def quadratic_polynomial_solution_2d_with_A_eps_eminus2():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]"
    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]"
    d2u_dy2_expr = "(-2)*(1 - x[0])*x[0]"


    # Coefficients
    eps = 1e-2
    A_expr = [["eps",  "0"], ["0", "eps"]]
    min_eig_A = eps
    a_expr = ["1.0", "1.0"]
    #lambda_expr = "1.0"
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

# checked
def example_1_sicom_paper():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]"
    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]"
    d2u_dy2_expr = "(-2)*(1 - x[0])*x[0]"


    # Coefficients
    eps = 0.1
    A_expr = [["eps",  "0"], ["0", "eps"]]
    min_eig_A = eps
    a_expr = ["2.0", "1.0"]
    #lambda_expr = "1.0"
    lambda_expr = "1.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

# checked
def example_2_sicom_paper():

    # Solution and it's derivatives
    u_expr = "(x[0] + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - x[1])*x[1]"
    du_dx_expr = "(1 + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - x[1])*x[1]"
    du_dy_expr = "(x[0] + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - 2*x[1])"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "((exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - x[1])*x[1]"
    d2u_dy2_expr = "(x[0] + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(-2)"

    # Coefficients
    eps = 1e-3
    A_expr = [["eps", "0"], ["0", "eps"]]
    min_eig_A = 1e-3
    a_expr = ["1.0", "0.0"]
    #lambda_expr = "1.0"
    lambda_expr = "1.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


# not-checked
def example_3_sicom_paper():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - x[2])*x[2]"
    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]*(1 - x[2])*x[2]"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])*(1 - x[2])*x[2]"
    du_dz_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])*(1 - 2*x[2])"
    grad_u_expr = [du_dx_expr, du_dy_expr, du_dz_expr]
    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]*(1 - x[2])*x[2]"
    d2u_dy2_expr = "(-2)*(1 - x[0])*x[0]*(1 - x[2])*x[2]"
    d2u_dz2_expr = "(-2)*(1 - x[0])*x[0]*(1 - x[1])*x[1]"


    # Coefficients
    eps = 0.01
    A_expr = [[eps,  0, 0], [0, eps, 0], [0, 0, eps]]
    min_eig_A = eps
    a_expr = ["3.0", "2.0", "1.0"]
    #lambda_expr = "1.0"
    lambda_expr = "10.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")" + "+" + \
                      "(" + d2u_dz2_expr + ")" + "*" + "(" + str(A_expr[2][2]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")" + "+" \
                 + "(" + a_expr[2] + ")" + "*" + "(" + grad_u_expr[2] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 3

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


# not-checked
def example_4_sicom_paper():

    # Solution and it's derivatives
    u_expr = "(x[0] + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - x[1])*x[1]*(1 - x[2])*x[2]"
    du_dx_expr = "(1 + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - x[1])*x[1]*(1 - x[2])*x[2]"
    du_dy_expr = "(x[0] + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - 2*x[1])*(1 - x[2])*x[2]"
    du_dz_expr = "(x[0] + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - x[1])*x[1]*(1 - 2*x[2])"
    grad_u_expr = [du_dx_expr, du_dy_expr, du_dz_expr]
    d2u_dx2_expr = "((exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - x[1])*x[1]"
    d2u_dy2_expr = "(x[0] + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(-2)"
    d2u_dz2_expr = "(x[0] + (exp((x[0] - 1)/1e-3) - exp(-1/1e-3))/(exp(-1/1e-3) - 1))*(1 - x[1])*x[1]*(- 2)"
    
    # Coefficients
    eps = 0.001
    A_expr = [["eps",  "0", "0"], ["0", "eps", "0"], ["0", "0", "eps"]]
    min_eig_A = eps
    a_expr = ["3.0", "2.0", "1.0"]
    #lambda_expr = "1.0"
    lambda_expr = "10.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")" + "+" + \
                      "(" + d2u_dz2_expr + ")" + "*" + "(" + str(A_expr[2][2]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")" + "+" \
                 + "(" + a_expr[2] + ")" + "*" + "(" + grad_u_expr[2] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 3

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

def example_1_iga_majorant_test():

    f_expr = "(pi*pi + 4*pi*pi)*sin(pi*x[0])*sin(2*pi*x[1])"
    eps = 1
    A_expr = [["eps", "0"], ["0", "eps"]]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["0", "0"]

    # Dirichlet BC
    uD_expr = "sin(pi*x[0])*sin(2*pi*x[1])"
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = "sin(pi*x[0])*sin(2*pi*x[1])"
    du_dx_expr = "pi*cos(pi*x[0])*sin(2*pi*x[1])"
    du_dy_expr = "2*pi*sin(pi*x[0])*cos(2*pi*x[1])"
    grad_u_expr = [du_dx_expr, du_dy_expr]


    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def long_chen_example_1_2d():

    eps = 1e-8
    
    # Solution and it's derivatives
    u_expr = "(x[0]*x[0] - exp((x[0] - 1)/eps))*x[1]*(1 - x[1])"
    du_dx_expr = "(2*x[0] - 1/eps*exp((x[0] - 1)/eps) )*x[1]*(1 - x[1])"
    du_dy_expr = "(x[0]*x[0] - exp((x[0] - 1)/eps) )*(1 - 2*x[1])"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "(2 - 1/(eps*eps)*exp((x[0] - 1)/eps) )*x[1]*(1 - x[1])"
    d2u_dy2_expr = "(x[0]*x[0] - exp((x[0] - 1)/eps) )*(-2)"

    # Coefficients
    A_expr = [["eps", "0"], ["0", "eps"]]
    A_inv_expr = [["1/eps", "0"], ["0", "1/eps"]]
    min_eig_A = eps
    a_expr = ["1.0", "0.0"]
    lambda_expr = "0.0"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + A_expr[0][0] + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + A_expr[1][1] + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, A_inv_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def long_chen_example_2_2d():

    eps = 1e-8
    
    # Solution and it's derivatives
    u_expr = "x[0]*x[0]*( x[1]*(1-x[1])+exp(-x[1]/sqrt(eps)))+exp((x[1]-1)/sqrt(eps))) )"
    du_dx_expr = "2*x[0]*( x[1]*(1-x[1])+exp(-x[1]/sqrt(eps)))+exp((x[1]-1)/sqrt(eps))) )"
    du_dy_expr = "x[0]*x[0]*( 1-2*x[1]-1/sqrt(eps))*exp(-x[1]/sqrt(eps)))+1/sqrt(eps))*exp((x[1]-1)/sqrt(eps))) )"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "2*( x[1]*(1-x[1])+exp(-x[1]/sqrt(eps)))+exp((x[1]-1)/sqrt(eps))) )"
    d2u_dy2_expr = "x[0]*x[0]*( -2+1/eps)*exp(-x[1]/sqrt(eps)))+1/eps)*exp((x[1]-1)/sqrt(eps))) )"

    # Coefficients
    A_expr = [["eps", "0"], ["0", "eps"]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1]])
    a_expr = ["1.0", "0.0"]
    lambda_expr = "0.0"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"
    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag


def long_chen_example_3_2d():

    eps = 1e-4
    # Solution and it's derivatives
    u_expr = "x[0]*x[1]*x[1] - x[1]*x[1]*exp(2*(x[0]-1)/eps)) - x[0]*exp(3*(x[1]-1)/eps)) " \
             "+ exp((2*(x[0]-1)+3*(x[1]-1))/eps))"
    du_dx_expr = "x[1]*x[1] - 2/eps)*x[1]*x[1]*exp(2*(x[0]-1)/eps)) - exp(3*(x[1]-1)/eps)) + 2/eps)*exp((2*(x[0]-1)+3*(x[1]-1))/eps))"
    du_dy_expr = "2*x[0]*x[1] - 2*x[1]*exp(2*(x[0]-1)/eps)) - 3/eps)*x[0]*exp(3*(x[1]-1)/eps)) + 3/eps)*exp((2*(x[0]-1)+3*(x[1]-1))/eps))"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "- 2*2/(" + str(eps*eps) + ") * x[1]*x[1]*exp(2*(x[0]-1)/eps)) + (2*2/(" + str(eps*eps) + ")) * exp((2*(x[0]-1)+3*(x[1]-1))/eps))"
    d2u_dy2_expr = "2*x[0] - 2*exp(2*(x[0]-1)/eps)) - (3*3/(" + str(eps*eps) + ")) * x[0]*exp(3*(x[1]-1)/eps)) + (3*3/(" + str(eps*eps) + ")) *exp((2*(x[0]-1)+3*(x[1]-1))/eps))"

    # Coefficients
    A_expr = [["eps", "0"], ["0", "eps"]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1]])
    a_expr = ["2.0", "3.0"]
    lambda_expr = "1.0"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    #divA_gradu_expr = "0.0"
    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"

    uD_expr = u_expr
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "reaction-convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag




def example_bank():

    f_expr = "x[0]"
    eps = 1.0
    A_expr = [["eps", "0"], ["0", "eps"]]
    min_eig_A = eps
    lambda_expr = "0.0"
    a_expr = ["cos(2*pi*x[1])", "1"]

    # Dirichlet BC
    uD_expr = "0"
    uN_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"
    pde_tag = "convection-diffusion-pde"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, eps, \
           uD_expr, uN_expr, domain, dim, solution_tag, material_tag, pde_tag

