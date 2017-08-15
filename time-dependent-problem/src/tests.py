__author__ = 'svetlana'


#-----------------------------------------------------------------------#
# test: unit domain in R^1, solution polynomial in space and time
#-----------------------------------------------------------------------#
# checked
def polynomial_in_space_time_solution_1d_t():

    # Solution and it's derivatives
    u_expr = '(1 - x[0])*x[0]*(t*t + t + 1)'
    du_dx_expr = '(1 - 2 * x[0])*(t*t + t + 1)'
    u_t_expr = '(1 - x[0])*x[0]*(2*t + 1)'
    grad_u_expr = [du_dx_expr]
    d2u_dx2_expr = '(- 2)*(t*t + t + 1)'

    # Coefficients
    A_expr = "1.0"
    min_eig_A = 1.0
    a_expr = "0.0"
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + A_expr + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"

    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 1.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

#-----------------------------------------------------------------------#
# test: unit domain in R^1, solution polynomial in space and
# trigonometric in time
#-----------------------------------------------------------------------#
# checked
def polynomial_in_space_trigonometric_in_time_solution_1d_t():

    # Solution and it's derivatives
    u_expr = '(1 - x[0])*x[0]*(sin(t) + cos(t))'
    du_dx_expr = '(1 - 2 * x[0])*(sin(t) + cos(t))'
    u_t_expr = '(1 - x[0])*x[0]*(cos(t) - sin(t))'
    d2u_dx2_expr = '(-2)*(sin(t) + cos(t))'
    grad_u_expr = [du_dx_expr]

    # Coefficients
    A_expr = "1.0"
    min_eig_A = 1.0
    a_expr = "0.0"
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + A_expr + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"

    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 1.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

#-----------------------------------------------------------------------#
# test: unit domain in R^1, solution polynomial in space and exponential
# in time
#-----------------------------------------------------------------------#
#
def polynomial_in_space_exponential_in_time_1d_t():

    # Solution and it's derivatives
    u_expr = '(1 - x[0])*x[0]*exp(t)'
    du_dx_expr = '(1 - 2 * x[0])*exp(t)'
    u_t_expr = '(1 - x[0])*x[0]*exp(t)'
    d2u_dx2_expr = '(-2)*exp(t)'
    grad_u_expr = [du_dx_expr]

    # Coefficients
    A_expr = "1.0"
    min_eig_A = 1.0
    a_expr = "0.0"
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + A_expr + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"

    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 1.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

#-----------------------------------------------------------------------#
# test: unit domain in R^1, solution polynomial in space and linear
# in time
#-----------------------------------------------------------------------#
# checked
def polynomial_in_space_linear_in_time_1d_t():

    # Solution and it's derivatives
    u_expr = '(1 - x[0])*x[0]*(t + 1)'
    du_dx_expr = '(1 - 2 * x[0])*(t + 1)'
    u_t_expr = '(1 - x[0])*x[0]'
    d2u_dx2_expr = '(-2)*(t + 1)'
    grad_u_expr = [du_dx_expr]

    # Coefficients
    A_expr = "1.0"
    min_eig_A = 1.0
    a_expr = "0.0"
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + A_expr + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"

    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 1.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag


#-----------------------------------------------------------------------#
# test: l-shape domain  in R^2, solution with singularities
#-----------------------------------------------------------------------#
# checked
def solution_with_singularities_2d_t():
    
    # Solution and it's derivatives
    u_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "pow(pow(x[0], 2) + pow(x[1], 2), 1.0 / 3.0) * " \
             "sin(2.0 / 3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))) * " \
             "(t*t + t + 1)"
    grad_u_x_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    " - x[1] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.)))) * " \
                    "(t*t + t + 1)"
    grad_u_y_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    " + x[1] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.)))) * " \
                    "(t*t + t + 1)"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]

    # Coefficients
    A_expr = "1.0"
    min_eig_A = 1.0
    a_expr = "0.0"
    lambda_expr = "0.0"

    f_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "pow(pow(x[0], 2) + pow(x[1], 2), 1.0/3.0) * " \
             "sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) * (2*t + 1)"

    domain = "l-shape-domain"
    dim = 2
    T = 1.0

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    uD_expr = u_expr

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

#-----------------------------------------------------------------------#
# test: l-shape domain  in R^2, solution with singularities, static
#-----------------------------------------------------------------------#

def static_solution_with_singularities_2d_t():
    # Solutions and its derivatives
    u_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "pow(pow(x[0], 2) + pow(x[1], 2), 1.0 / 3.0) * " \
             "sin(2.0 / 3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.)))"
    grad_u_x_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    " - x[1] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))))"
    grad_u_y_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    " + x[1] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))))"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]

    # Coefficients
    A_expr = "1.0"
    min_eig_A = 1.0
    a_expr = "0.0"
    lambda_expr = "0.0"

    f_expr = "0.0"

    domain = "l-shape-domain"
    dim = 2
    T = 1.0

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    uD_expr = u_expr

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

#-----------------------------------------------------------------------#
# test: l-shape domain  in R^2, solution without singularities
#-----------------------------------------------------------------------#

def solution_without_singularities_2d_t():

    # Solution and it's derivatives
    u_expr = "cos(x[0])*exp(x[1])*(t*t + t + 1)"
    grad_u_x_expr = "-sin(x[0]) * exp(x[1])*(t*t + t + 1)"
    grad_u_y_expr = "cos(x[0]) * exp(x[1])*(t*t + t + 1)"
    u_t_expr = 'cos(x[0])*exp(x[1])*(2*t + 1)'
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]
    d2u_dx2_expr = "0.0"
    d2u_dy2_expr = "cos(x[0])*exp(x[1])*(t*t + t + 1)"

    # Coefficients
    A_expr = [[1,  0], [0, 1]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1]])
    a_expr = ["0.0", "0.0"]
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"


    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
                   + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                   + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                   + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"

    uD_expr = u_expr

    domain = "pi-shape-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 1.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

#-----------------------------------------------------------------------#
# test: unit domain in R^3, polynomial solution
#-----------------------------------------------------------------------#
# checked
def polynomial_solution_3d_t():

    u_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - x[2])*x[2]*(t*t + t + 1)'

    grad_u_x_expr = '(1 - 2*x[0])*(1 - x[1])*x[1]*(1 - x[2])*x[2]*(t*t + t + 1)'
    grad_u_y_expr = '(1 - x[0])*x[0]*(1 - 2*x[1])*(1 - x[2])*x[2]*(t*t + t + 1)'
    grad_u_z_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - 2*x[2])*(t*t + t + 1)'
    u_t_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - x[2])*x[2]*(2*t + 1)'
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr, grad_u_z_expr]
    d2u_dx2_expr = '(- 2)*(1 - x[1])*x[1]*(1 - x[2])*x[2]*(t*t + t + 1)'
    d2u_dy2_expr = '(1 - x[0])*x[0]*(- 2)*(1 - x[2])*x[2]*(t*t + t + 1)'
    d2u_dz2_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(- 2)*(t*t + t + 1)'

    # Coefficients
    A_expr = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1], A_expr[2][2]])
    a_expr = ["0.0", "0.0", "0.0"]
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")" + "+" + \
                      "(" + d2u_dz2_expr + ")" + "*" + "(" + str(A_expr[2][2]) + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
             + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")" + "+" \
             + "(" + a_expr[2] + ")" + "*" + "(" + grad_u_expr[2] + ")"

    uD_expr = u_expr

    domain = "unit-domain"
    dim = 3

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 10.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

#-----------------------------------------------------------------------#
# test: l-shape domain  in R^3, solution with singularities
#-----------------------------------------------------------------------#

def solution_with_singularities_3d_t():
    u_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "(pow(pow(x[0], 2) + pow(x[1], 2), 1.0 / 3.0) * " \
             "sin(2.0 / 3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))) + x[2]) * " \
             "(t*t + t + 1)"
    grad_u_x_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "(2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    "- x[1] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.)))) + x[2]) * (t*t + t + 1)"
    grad_u_y_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "(2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    "+ x[1] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.)))) + x[2]) * (t*t + t + 1)"
    grad_u_z_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "(t*t + t + 1)"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr, grad_u_z_expr]

    f_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
        "(pow(pow(x[0], 2) + pow(x[1], 2), 1.0/3.0) * " \
        "sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) + x[2]) * " \
        "(2*t + 1)"

    # Coefficients
    A_expr = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1], A_expr[2][2]])
    a_expr = ["0.0", "0.0", "0.0"]
    lambda_expr = "0.0"

    domain = "l-shape-domain"
    dim = 3
    T = 1.0

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

#-----------------------------------------------------------------------#
# test: l-shape domain in R^3, solution without singularities
#-----------------------------------------------------------------------#

def solution_without_singularities_3d_t():
    u_expr = "cos(x[0])*exp(x[1])*(1 - x[2])*(t*t + t + 1)"
    grad_u_x_expr = "-sin(x[0]) * exp(x[1])*(1 - x[2])*(t*t + t + 1)"
    grad_u_y_expr = "cos(x[0]) * exp(x[1])*(1 - x[2])*(t*t + t + 1)"
    grad_u_z_expr = "cos(x[0])*exp(x[1])*(-1)*(t*t + t + 1)"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr, grad_u_z_expr]
    u_t_expr = 'cos(x[0])*exp(x[1])*(1 - x[2])*(2*t + 1)'

    d2u_dx2_expr = '-cos(x[0]) * exp(x[1])*(1 - x[2])*(t*t + t + 1)'
    d2u_dy2_expr = 'cos(x[0]) * exp(x[1])*(1 - x[2])*(t*t + t + 1)'
    d2u_dz2_expr = '0.0'

    # Coefficients
    A_expr = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1], A_expr[2][2]])
    a_expr = ["0.0", "0.0", "0.0"]
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")" + "+" + \
                      "(" + d2u_dz2_expr + ")" + "*" + "(" + str(A_expr[2][2]) + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
             + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")" + "+" \
             + "(" + a_expr[2] + ")" + "*" + "(" + grad_u_expr[2] + ")"

    uD_expr = u_expr

    domain = "unit-domain"
    dim = 3

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 10.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag


#-----------------------------------------------------------------------#
# test: unit domain in R^2, polynomial solution
#-----------------------------------------------------------------------#

def polynomial_solution_2d_t():

    u_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(t*t + t + 1)'
    grad_u_x_expr = '(1 - 2 * x[0])*(1 - x[1])*x[1]*(t*t + t + 1)'
    grad_u_y_expr = '(1 - 2 * x[1])*(1 - x[0])*x[0]*(t*t + t + 1)'
    u_t_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(2*t + 1)'
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]
    d2u_dx2_expr = "(- 2)*(1 - x[1])*x[1]*(t*t + t + 1)"
    d2u_dy2_expr = "(- 2)*(1 - x[0])*x[0]*(t*t + t + 1)"

    # Coefficients
    A_expr = [[1.0, 0], [0, 1.0]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1]])
    a_expr = ["0.0", "0.0"]
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
             + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"

    uD_expr = u_expr

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 1.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag


#-----------------------------------------------------------------------#
# test: unit domain in R^2, trigonometric solution
# example 8 from the AMC paper
#-----------------------------------------------------------------------#

def trigonometric_2d_t():

    u_expr = '(sin(pi*x[0])*sin(3*pi*x[1]) + sin(pi*x[1])*sin(3*pi*x[0]))*(pow(t, 3) + sin(t) + 1)'
    grad_u_x_expr = '(pi*cos(pi*x[0])*sin(3*pi*x[1]) + 3*pi*cos(3*pi*x[0])*sin(pi*x[1]))*(pow(t, 3) + sin(t) + 1)'
    grad_u_y_expr = '(pi*cos(pi*x[1])*sin(3*pi*x[0]) + 3*pi*cos(3*pi*x[1])*sin(pi*x[0]))*(pow(t, 3) + sin(t) + 1)'
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]
    u_t_expr = '(sin(pi*x[0])*sin(3*pi*x[1]) + sin(pi*x[1])*sin(3*pi*x[0]))*(3*pow(t, 2) + cos(t))'
    d2u_dx2_expr = "(-pi*pi*sin(pi*x[0])*sin(3*pi*x[1]) - 9*pi*pi*sin(3*pi*x[0])*sin(pi*x[1]))*(pow(t, 3) + sin(t) + 1)"
    d2u_dy2_expr = "(-pi*pi*sin(pi*x[1])*sin(3*pi*x[0]) - 9*pi*pi*sin(3*pi*x[1])*sin(pi*x[0]))*(pow(t, 3) + sin(t) + 1)"

    # Coefficients
    A_expr = [[1, 0], [0, 1]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1]])
    a_expr = ["0.0", "0.0"]
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
             + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"

    uD_expr = u_expr

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 0.1

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

#-----------------------------------------------------------------------#
# test: unit domain in R^2, trigonometric and drasticly changing
# in time solution
# example 14 from the paper
#-----------------------------------------------------------------------#

def trigonometric_drastic_changing_in_time_2d_t():

    u_expr = 'sin(pi*x[0])*sin(3*pi*x[1])*(t*sin(t) + 1) ' \
             '+ sin(3*pi*x[0])*sin(pi*x[1])*(cos(t) + t*sin(3*t))'
    grad_u_x_expr = '3*pi*cos(3*pi*x[0])*sin(pi*x[1])*(cos(t) + t*sin(3*t)) ' \
                    '+ pi*cos(pi*x[0])*sin(3*pi*x[1])*(t*sin(t) + 1)'
    grad_u_y_expr = 'pi*cos(pi*x[1])*sin(3*pi*x[0])*(cos(t) + t*sin(3*t)) + 3*pi*cos(3*pi*x[1])*sin(pi*x[0])*(t*sin(t) + 1)'
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]

    u_t_expr = '(sin(pi*x[0])*sin(3*pi*x[1])*(sin(t) + t*cos(t)) + sin(3*pi*x[0])*sin(pi*x[1])*(sin(3*t) - sin(t) + 3*t*cos(3*t)) )'
    d2u_dx2_expr = "-10*pow(pi, 2)*sin(pi*x[0])*sin(3*pi*x[1])*(t*sin(t) + 1)"
    d2u_dy2_expr = "-10*pow(pi, 2)*sin(3*pi*x[0])*sin(pi*x[1])*(cos(t) + t*sin(3*t))"

    # Coefficients
    A_expr = [[1, 0], [0, 1]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1]])
    a_expr = ["0.0", "0.0"]
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
             + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"

    uD_expr = u_expr

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 0.1

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag


def solution_without_singularities_on_circle_2d_t():
    u_expr = "cos(x[0])*exp(x[1])*(t*t + t + 1)"
    grad_u_x_expr = "-sin(x[0]) * exp(x[1])*(t*t + t + 1)"
    grad_u_y_expr = "cos(x[0]) * exp(x[1])*(t*t + t + 1)"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]
    u_t_expr = "cos(x[0])*exp(x[1])*(2*t + 1)"

    d2u_dx2_expr = "-cos(x[0]) * exp(x[1])*(t*t + t + 1)"
    d2u_dy2_expr = "cos(x[0]) * exp(x[1])*(t*t + t + 1)"

    # Coefficients
    A_expr = [[1, 0], [0, 1]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1]])
    a_expr = ["0.0", "0.0"]
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
             + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"

    uD_expr = u_expr

    domain = "circle-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 0.1

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag


def solution_with_singularities_exp_unit_circle_2d_t():
    u_expr = "exp(x[0] + x[1])*(t*sin(t) + 1)"
    grad_u_x_expr = "exp(x[0] + x[1])*(t*sin(t) + 1)"
    grad_u_y_expr = "exp(x[0] + x[1])*(t*sin(t) + 1)"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]
    u_t_expr = "exp(x[0] + x[1])*(sin(t) + t*cos(t))"

    d2u_dx2_expr = "exp(x[0] + x[1])*(t*sin(t) + 1)"
    d2u_dy2_expr = "exp(x[0] + x[1])*(t*sin(t) + 1)"

    # Coefficients
    A_expr = [[1, 0], [0, 1]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1]])
    a_expr = ["0.0", "0.0"]
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
             + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"

    uD_expr = u_expr

    domain = "circle-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 3.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

# new example

def solution_oscilation_exp_unit_square_2d_t():
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*sin(17*x[0]*x[1])*exp(x[0] + x[1])*(t*sin(t) + 1)"
    grad_u_x_expr = "-sin(x[0]) * exp(x[1])*(t*t + t + 1)"
    grad_u_y_expr = "cos(x[0]) * exp(x[1])*(t*t + t + 1)"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]
    f_expr = "cos(x[0])*exp(x[1])*(2*t + 1)"

    domain = "unit-square-domain"
    dim = 2
    T = 3.0

    return u_expr, grad_u_expr, f_expr, domain, dim, T

    u_t_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*sin(17*x[0]*x[1])*exp(x[0] + x[1])*(sin(t) + t*cos(t))"

    d2u_dx2_expr = "exp(x[0] + x[1])*(t*sin(t) + 1)"
    d2u_dy2_expr = "exp(x[0] + x[1])*(t*sin(t) + 1)"

    # Coefficients
    A_expr = [[1, 0], [0, 1]]
    min_eig_A = min([A_expr[0][0], A_expr[1][1]])
    a_expr = ["0.0", "0.0"]
    lambda_expr = "0.0"

    # A * grad u = [(du_dx_expr * A_11 + du_dy_expr * A_12), (du_dx_expr * A_21 + du_dy_expr * A_22)]
    # div (A * grad u) = d2u_dx2_expr * A_11 + d2u_dy2_expr * A_22
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = u_t_expr \
             + "-" + "(" + divA_gradu_expr + ")" + "+" \
             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
             + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
             + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"

    uD_expr = u_expr

    domain = "circle-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    T = 3.0

    return u_expr, grad_u_expr, f_expr, A_expr, min_eig_A, lambda_expr, a_expr, uD_expr, domain, T, dim, solution_tag, material_tag

def solution_exp_singularity_unit_square_2d_t():
    u_expr = "(1 - x[0])*x[0]*" \
             "(1 - x[1])*x[1]*" \
             "exp(-10*(" \
             "     (x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.6)*(x[1] - 0.6)" \
             "          ))*(t*sin(t) + 1)"
    grad_u_x_expr = "-sin(x[0]) * exp(x[1])*(t*sin(t) + 1)"
    grad_u_y_expr = "cos(x[0]) * exp(x[1])*(t*sin(t) + 1)"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]
    f_expr = "cos(x[0])*exp(x[1])*(...)"

    domain = "unit-square-domain"
    dim = 2
    T = 3.0

    return u_expr, grad_u_expr, f_expr, domain, dim, T

def solution_exp_singularity_unit_square_2d_t():
    u_expr = "exp(-10*((x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.6)*(x[1] - 0.6)))*(t*sin(t) + 1)"
    grad_u_x_expr = "-20*(x[0] - 0.2)*exp(-10*((x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.6)*(x[1] - 0.6)))*(t*sin(t) + 1)"
    grad_u_y_expr = "-20*(x[1] - 0.6)*exp(-10*((x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.6)*(x[1] - 0.6)))*(t*sin(t) + 1)"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]
    f_expr = "40*exp(-10*((x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.6)*(x[1] - 0.6)))*(t*sin(t) + 1) + exp(-10*((x[0] - 0.2)*(x[0] - 0.2) + (x[1] - 0.6)*(x[1] - 0.6)))*(t*cos(t) + sin(t))"

    domain = "unit-square-domain"
    dim = 2
    T = 3.0

    return u_expr, grad_u_expr, f_expr, domain, dim, T

def solution_with_singularities_sint_t_2d_t():

    u_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "(pow(pow(x[0], 2) + pow(x[1], 2), 1.0 / 3.0) * " \
             "sin(2.0 / 3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.)))) * " \
             "(sin(t)*t + 1)"
    du_dx_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "(2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    "- x[1] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))))) * (sin(t)*t + 1)"
    du_dy_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "(2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    "+ x[1] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))))) * (sin(t)*t + 1)"
    du_dt_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "(2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    "+ x[1] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))))) * (sin(t) + cos(t)*t)"
    grad_u_expr = [du_dx_expr, du_dy_expr]

    f_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
        "(pow(pow(x[0], 2) + pow(x[1], 2), 1.0/3.0) * " \
        "sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.)))) * " \
        "(sin(t) + cos(t)*t)"

    domain = "l-shape-domain"
    dim = 2
    T = 1.0

    return u_expr, grad_u_expr, f_expr, domain, dim, T