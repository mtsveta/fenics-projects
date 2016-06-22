__author__ = "Svetlana Matculevich <svetlana.v.matculevich@jyu.fi>"
__date__ = ""
__copyright__ = ""
__license__ = ""

# checked
def quadratic_polynomial_solution_1d_sigma_1000():

    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["1"]
    lambda_1 = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 / (1000.0 * sqrt(2 * pi)) * exp( - pow(x[0] - 0.5, 2) / (2.0 * pow(1000.0, 2)))"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_1d_sigma_1():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["1"]
    lambda_1 = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 / (1.0 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(1.0, 2)))"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_1d_sigma_5minus2():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["1"]
    lambda_1 = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 / (0.05 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(0.05, 2)))"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def linear_trigonometric_solution_1d_rho_1minus3():

    # Solution and it's derivatives
    u_expr = "x[0]*sin(3*pi*x[0])"
    du_dx_expr = "sin(3*pi*x[0]) + x[0]*3*pi*cos(3*pi*x[0])"
    d2u_dx2_expr = "6*pi*cos(3*pi*x[0]) - 9*pow(pi, 2)*x[0]*sin(3*pi*x[0])"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["1"]
    lambda_1 = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "0.001 * (x[0] + 0.001)"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def linear_trigonometric_solution_1d_rho_1():

    # Solution and it's derivatives
    u_expr = "x[0]*sin(3*pi*x[0])"
    du_dx_expr = "sin(3*pi*x[0]) + x[0]*3*pi*cos(3*pi*x[0])"
    d2u_dx2_expr = "6*pi*cos(3*pi*x[0]) - 9*pow(pi, 2)*x[0]*sin(3*pi*x[0])"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["1"]
    lambda_1 = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 * (x[0] + 0.001)"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def linear_trigonometric_solution_1d_rho_1000():

    # Solution and it's derivatives
    u_expr = "x[0]*sin(3*pi*x[0])"
    du_dx_expr = "sin(3*pi*x[0]) + x[0]*3*pi*cos(3*pi*x[0])"
    d2u_dx2_expr = "6*pi*cos(3*pi*x[0]) - 9*pow(pi, 2)*x[0]*sin(3*pi*x[0])"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["1"]
    lambda_1 = 1.0
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1000.0 * (x[0] + 0.001)"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked

def quadratic_polynomial_solution_2d_sigma_x_1_sigma_y_5():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]"
    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]"
    d2u_dy2_expr = "(-2)*(1 - x[0])*x[0]"

    # Coefficients
    A_expr = [[1,  0], [0, 1]]
    lambda_1 = 1.0
    a_expr = ["0", "0"]
    lambda_expr = "1.0 / (1.0 * 5.0 * sqrt(2 * pi)) * exp( - 0.5 * ( pow(x[0] - 0.5, 2) / pow(1.0, 2) + " \
                                                                    "pow(x[1] - 0.5, 2) / pow(5.0, 2) - " \
                                                                    "(x[1] - 0.5) * (x[1] - 0.5) / (1.0 * 5.0) ) )"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_2d_sigma_x_1minus2_sigma_y_5():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]"
    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]"
    d2u_dy2_expr = "(-2)*(1 - x[0])*x[0]"


    # Coefficients
    A_expr = [[1,  0], [0, 1]]
    lambda_1 = 1.0
    a_expr = ["0", "0"]
    lambda_expr = "1.0 / (0.01 * 5.0 * sqrt(2 * pi)) * exp( - 0.5 * ( pow(x[0] - 0.5, 2) / pow(1.0, 2) + " \
                                                                    "pow(x[1] - 0.5, 2) / pow(5.0, 2) - " \
                                                                    "(x[1] - 0.5) * (x[1] - 0.5) / (0.01 * 5.0) ) )"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def linear_trigonometric_solution_2d_rho_1000():

    # Solution and it's derivatives
    u_expr = "x[0]*sin(3*pi*x[0])*x[1]*sin(pi*x[1])"
    du_dx_expr = "(sin(3*pi*x[0]) + x[0]*3*pi*cos(3*pi*x[0]))*x[1]*sin(pi*x[1])"
    du_dy_expr = "(sin(pi*x[1]) + x[1]*pi*cos(pi*x[1]))*x[0]*sin(3*pi*x[0])"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "(6*pi*cos(3*pi*x[0]) - 9*pow(pi, 2)*x[0]*sin(3*pi*x[0]))*x[1]*sin(pi*x[1])"
    d2u_dy2_expr = "(2*pi*cos(pi*x[1]) - x[1]*pi*pi*sin(pi*x[1]))*x[0]*sin(3*pi*x[0])"

    # Coefficients
    A_expr = [[1,  0], [0, 1]]
    lambda_1 = 1.0
    a_expr = ["0", "0"]
    lambda_expr = "1000.0 * (x[0] + 0.01) * (x[1] + 0.01)"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"


    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_2d_with_A():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]"
    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]"
    d2u_dy2_expr = "(-2)*(1 - x[0])*x[0]"


    # Coefficients
    A_expr = [[1,  0], [0, 10]]
    lambda_1 = 1.0
    a_expr = ["-1.0", "1.0"]
    #lambda_expr = "1.0"
    lambda_expr = "1.0 / (1.0 * 5.0 * sqrt(2 * pi)) * exp( - 0.5 * ( pow(x[0] - 0.8, 2) / pow(1.0, 2) + " \
                                                                    "pow(x[1] - 0.3, 2) / pow(5.0, 2) - " \
                                                                    "(x[1] - 0.8) * (x[1] - 0.3) / (1.0 * 5.0) ) )"

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

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_2d_with_A_sigma_5minus2():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]"
    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])"
    grad_u_expr = [du_dx_expr, du_dy_expr]
    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]"
    d2u_dy2_expr = "(-2)*(1 - x[0])*x[0]"

    # Coefficients
    A_expr = [[1,  0], [0, 10]]
    lambda_1 = 1.0
    a_expr = ["-1.0", "1.0"]
    lambda_expr = "1.0 / (0.05 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(0.05, 2)))"


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

    domain = "unit-domain"
    dim = 2

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_3d_with_A():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - x[2])*x[2]"
    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]*(1 - x[2])*x[2]"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])*(1 - x[2])*x[2]"
    du_dz_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - 2*x[2])"
    grad_u_expr = [du_dx_expr, du_dy_expr, du_dz_expr]
    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]*(1 - x[2])*x[2]"
    d2u_dy2_expr = "(1 - x[0])*x[0]*(-2)*(1 - x[2])*x[2]"
    d2u_dz2_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*(-2)"

    A_expr = [[1, 0, 0], [0, 2, 0], [0, 0, 5]]
    lambda_1 = 1.0
    lambda_expr = "0.0"
    a_expr = ["-x[0] + 1", "1", "1"]

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

    domain = "unit-domain"
    dim = 3

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_1d_sigma_1_epsilon_1eminus2():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["0.01"]
    lambda_1 = 0.01
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 / (1.0 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(1.0, 2)))"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_1d_sigma_1minus2():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["1"]
    lambda_1 = 1
    a_expr = ["-x[0] + 1"]
    lambda_expr = "1.0 / (0.01 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(0.01, 2)))"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_1d_lmdb_zero():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["1"]
    lambda_1 = 1
    a_expr = ["-1"]
    #lambda_expr = "1.0 / (1.0 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(1.0, 2)))"
    lambda_expr = "0.0"
    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u
    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_1d_lmdb_zero_a_linear():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["10"]
    lambda_1 = 10
    a_expr = ["-x[0]"]
    lambda_expr = "0.0"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0]) + ")"

    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def quadratic_polynomial_solution_1d_lmdb_zero_a_cubic():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients

    A_expr = ["1"]
    lambda_1 = 1
    a_expr = ["-100*x[0]*x[0]*x[0]"]
    lambda_expr = "0.0"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0]) + ")"

    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_1d_lmdb_zero_a_cubic_epsilon_1minus2():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["0.01"]
    lambda_1 = 0.01
    a_expr = ["-x[0]*x[0]*x[0]"]
    lambda_expr = "0.0"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0]) + ")"

    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr
    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def quadratic_polynomial_solution_1d_lmdb_zero_eps_1minus2():

    # Solution and it's derivatives
    u_expr = "(1 - x[0])*x[0]"
    du_dx_expr = "(1 - 2*x[0])"
    d2u_dx2_expr = "-2"
    grad_u_expr = [du_dx_expr]

    # Coeffitients
    A_expr = ["1e-2"]
    lambda_1 = 1e-2
    a_expr = ["-1"]
    lambda_expr = "0.0"

    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0]) + ")"

    f_expr = "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def classic_example_1d_b_1_eps_1():

    eps = "1.0"
    b = "1.0"
    u_expr = "(1 - exp((" + b + ") * x[0] / (" + eps + ")))/(1 - exp(" + "(" + b + ") / (" + eps + ")))"
    du_dx_expr = "-(" + b + ") / (" + eps + ") * exp((" + b + ") * x[0] / (" + eps + "))/( 1 - exp(" + "(" + b + ") / (" + eps + ")))"
    grad_u_expr = [du_dx_expr]
    d2u_dx2_expr = "- (" + b + ") / (" + eps + ") * " + "(" + b + ") / (" + eps + ") * exp((" + b + ") * x[0] / (" + eps + "))/( 1 - exp(" + "(" + b + ") / (" + eps + ")))"

    A_expr = [eps]
    lambda_1 = 1.0
    lambda_expr = "0.0"
    a_expr = [b]

    f_expr = "0.0"
    uD_expr = u_expr
    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def classic_example_1d_b_1_eps_1eminis2():

    eps = "0.01"
    b = "1.0"
    u_expr = "(1 - exp((" + b + ") * x[0] / (" + eps + ")))/(1 - exp(" + "(" + b + ") / (" + eps + ")))"
    du_dx_expr = "-(" + b + ") / (" + eps + ") * exp((" + b + ") * x[0] / (" + eps + "))/( 1 - exp(" + "(" + b + ") / (" + eps + ")))"
    grad_u_expr = [du_dx_expr]
    d2u_dx2_expr = "- (" + b + ") / (" + eps + ") * " + "(" + b + ") / (" + eps + ") * exp((" + b + ") * x[0] / (" + eps + "))/( 1 - exp(" + "(" + b + ") / (" + eps + ")))"

    A_expr = [eps]
    lambda_1 = 0.01
    lambda_expr = "0.0"
    a_expr = [b]

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "0.0"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

# checked
def classic_example_1d_b_3_eps_1eminis2():

    eps = "0.001"
    b = "1.0"
    u_expr = "(1 - exp((" + b + ") * x[0] / (" + eps + ")))/(1 - exp(" + "(" + b + ") / (" + eps + ")))"
    du_dx_expr = "-(" + b + ") / (" + eps + ") * exp((" + b + ") * x[0] / (" + eps + "))/( 1 - exp(" + "(" + b + ") / (" + eps + ")))"
    grad_u_expr = [du_dx_expr]
    d2u_dx2_expr = "- (" + b + ") / (" + eps + ") * " + "(" + b + ") / (" + eps + ") * exp((" + b + ") * x[0] / (" + eps + "))/( 1 - exp(" + "(" + b + ") / (" + eps + ")))"

    A_expr = [eps]
    lambda_1 = 0.001
    lambda_expr = "0.0"
    a_expr = [b]

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "0.0"
    uD_expr = u_expr

    domain = "unit-domain"
    dim = 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag


def f_1_uD_polynomial_1d_sigma_1000():

    # Input data
    f_expr = "1.0"
    A_expr = ["1.0"]
    lambda_1 = 1.0
    lambda_expr = "1.0 / (1000.0 * sqrt(2 * pi)) * exp( - pow(x[0] - 0.5, 2) / (2.0 * pow(1000.0, 2)))"
    a_expr = ["0.0"]

    # Dirichlet BC
    uD_expr = "(1 - x[0])*x[0]"

    # Domain
    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_polynomial_1d_sigma_1():

    # Input data
    f_expr = "1.0"
    A_expr = ["1.0"]
    lambda_1 = 1.0
    lambda_expr = "1.0 / (1.0 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(1.0, 2)))"
    a_expr = ["0.0"]

    # Dirichlet BC
    uD_expr = "(1 - x[0])*x[0]"

    # Domain
    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_polynomial_1d_sigma_1minus2():

    # Input data
    f_expr = "1.0"
    A_expr = ["1.0"]
    lambda_1 = 1.0
    lambda_expr = "1.0 / (0.01 * sqrt(2 * pi)) * exp( -pow(x[0] - 0.5, 2) / (2.0 * pow(0.01, 2)))"
    a_expr = ["0.0"]

    # Dirichlet BC
    uD_expr = "(1 - x[0])*x[0]"

    # Domain
    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_polynomial_1d_lmdb_zero():

    # Input data
    f_expr = "1.0"
    A_expr = ["1.0"]
    lambda_1 = 1.0
    lambda_expr = "0.0"
    a_expr = ["0.0"]

    # Dirichlet BC
    uD_expr = "(1 - x[0])*x[0]"

    # Domain
    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_zero_1d_lmdb_zero_a_const():

    # Input data
    f_expr = "1.0"
    A_expr = ["0.01"]
    lambda_1 = 0.01
    lambda_expr = "0.0"
    a_expr = ["-1.0"]

    # Dirichlet BC
    uD_expr = "0.0"

    # Domain
    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_polynomial_1d_lmdb_zero_a_const_eps_1minus3():

    # Input data
    f_expr = "1.0"
    A_expr = ["1e-3"]
    lambda_1 = 1e-3
    lambda_expr = "0.0"
    a_expr = ["1.0"]

    # Dirichlet BC
    uD_expr = "(1 - x[0])*x[0]"

    # Domain
    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_zero_1d_lmdb_zero_a_cubic_eps_1minus3():

    # Input data
    f_expr = "1.0"
    A_expr = ["1e-3"]
    lambda_1 = 1e-3
    lambda_expr = "0.0"
    a_expr = ["-(x[0] + 1.0)*(x[0] + 1.0)*(x[0] + 1.0)"]

    # Dirichlet BC
    uD_expr = "(1 - x[0])*x[0]"

    # Domain
    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_zero_2d_lmdb_zero_a_const_eps_1minus2():

    # Input data
    f_expr = "1.0"
    eps = 1e-2
    A_expr = [[eps,  0], [0, eps]]
    lambda_1 = eps
    lambda_expr = "0.0"
    a_expr = ["1", "0"]

    # Dirichlet BC
    uD_expr = "0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_zero_2d_lmdb_zero_a_const_eps_1minus3():

    f_expr = "1.0"
    eps = 1e-3
    A_expr = [[eps,  0], [0, eps]]
    lambda_1 = eps
    lambda_expr = "0.0"
    a_expr = ["1", "1"]

    # Dirichlet BC
    uD_expr = "0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_trigonometric_uD_0_1d_rho_1minus3():


    f_expr = "x[0]*sin(3*pi*x[0])"
    A_expr = ["1.0"]
    lambda_1 = 1.0
    lambda_expr = "0.001 * (x[0] + 0.001)"
    a_expr = ["-1.0"]

    # Dirichlet BC
    uD_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_trigonometric_uD_0_1d_rho_1():

    f_expr = "x[0]*sin(3*pi*x[0])"
    A_expr = ["1.0"]
    lambda_1 = 1.0
    lambda_expr = "1.0 * (x[0] + 0.001)"
    a_expr = ["0.0"]

    # Dirichlet BC
    uD_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_trigonometric_uD_0_1d_rho_1000():

    f_expr = "x[0]*sin(3*pi*x[0])"
    A_expr = ["1.0"]
    lambda_1 = 1.0
    lambda_expr = "1000.0 * (x[0] + 0.001)"
    a_expr = ["0.0"]

    # Dirichlet BC
    uD_expr = "0.0"

    domain = "unit-domain"
    dim = 1

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_0_2d_sigma_x_1_sigma_y_5():

    f_expr = "1.0"
    A_expr = [[1,  0], [0, 1]]
    lambda_1 = min([A_expr[0][0], A_expr[1][1]])
    lambda_expr = "1.0 / (1.0 * 5.0 * sqrt(2 * pi)) * exp( - 0.5 * ( pow(x[0] - 0.5, 2) / pow(1.0, 2) + " \
                  "pow(x[1] - 0.5, 2) / pow(5.0, 2) - " \
                  "(x[1] - 0.5) * (x[1] - 0.5) / (1.0 * 5.0) ) )"
    a_expr = ["0.0", "0.0"]

    uD_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_trigonometric_uD_0_2d_rho_1000():

    f_expr = "x[0]*sin(3*pi*x[0])"
    A_expr = [[1.0,  0], [0, 1.0]]
    lambda_1 = min([A_expr[0][0], A_expr[1][1]])
    lambda_expr = "1000.0 * (x[0] + 0.01) * (x[1] + 0.01)"
    a_expr = ["0.0", "0.0"]

    uD_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_0_2d_with_A():

    f_expr = "1.0"
    A_expr = [[1,  0], [0, 10]]
    lambda_1 = min([A_expr[0][0], A_expr[1][1]])
    lambda_expr = "0.0"
    a_expr = ["0.0", "0.0"]

    uD_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_0_2d_with_A_sigma_x_5minus2_sigma_y_5():

    f_expr = "1.0"
    A_expr = [[1,  0.3], [0.3, 10]]
    lambda_1 = 0.990011
    lambda_expr = "1.0 / (0.05 * 5.0 * sqrt(2 * pi)) * exp( - 0.5 * ( pow(x[0] - 0.5, 2) / pow(0.05, 2) + " \
                  "pow(x[1] - 0.5, 2) / pow(5.0, 2) - " \
                  "(x[0] - 0.5) * (x[1] - 0.5) / (0.05 * 5.0) ) )"
    a_expr = ["0.0", "0.0"]

    uD_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def f_1_uD_0_2d_with_changing_A_lmdb_0():

    f_expr = "1.0"
    A_expr_1 = [[1.0, 0.0], [0.0, 10.0]]
    A_expr_2 = [[5.0, 0.0], [0.0, 1.0]]
    A_expr = [A_expr_1, A_expr_2]
    lambda_1 = min([A_expr_1[0][0], A_expr_1[1][1], A_expr_2[0][0], A_expr_2[1][1]])
    lambda_expr = "0.0"
    a_expr = ["0.0", "0.0"]

    uD_expr = "0.0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def quadratic_polynomial_solution_3d_with_A():

    f_expr = "1.0"
    A_expr = [[1, 0, 0], [0, 2, 0], [0, 0, 5]]
    lambda_1 = min([A_expr[0][0], A_expr[1][1], A_expr[2][2]])
    lambda_expr = "0.0"
    a_expr = ["0.0", "0.0", "0.0"]

    uD_expr = "1.0"

    domain = "unit-domain"
    dim = 3

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def solution_with_singularities_2d():

    f_expr = "0.0"
    A_expr = [[1, 0], [0, 1]]
    lambda_1 = min([A_expr[0][0], A_expr[1][1]])
    lambda_expr = "0.0"
    a_expr = ["0.0", "0.0"]

    uD_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
              "pow(pow(x[0], 2) + pow(x[1], 2), 1.0 / 3.0) * " \
              "sin(2.0 / 3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.)))"

    domain = "l-shape-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def solution_with_singularities_3d():

    f_expr = "0.0"
    A_expr = [[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]]
    lambda_1 = min([A_expr[0][0], A_expr[1][1], A_expr[2][2]])
    lambda_expr = "0.0"
    a_expr = ["0.0", "0.0", "0.0"]


    uD_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "pow(pow(x[0], 2) + pow(x[1], 2), 1.0 / 3.0) * " \
             "sin(2.0 / 3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))) " \
             "+ x[2]"

    domain = "l-shape-domain"
    dim = 3

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def example_1_book_stynes_2d():

    f_expr = "1.0"
    eps = 1e-3
    A_expr = [[eps,  0], [0, eps]]
    lambda_1 = eps
    lambda_expr = "0.0"
    a_expr = ["1", "0"]

    # Dirichlet BC
    uD_expr = "0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def example_2_eps_1eminus3book_stynes_2d():

    f_expr = "x[0]"
    eps = 1e-3
    A_expr = [[eps,  0], [0, eps]]
    lambda_1 = eps
    lambda_expr = "0.0"
    a_expr = ["-x[0]", "0"]

    # Dirichlet BC
    uD_expr = "0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def example_2_eps_1eminus2book_stynes_2d():

    f_expr = "x[0]"
    eps = 1e-2
    A_expr = [[eps,  0], [0, eps]]
    lambda_1 = eps
    lambda_expr = "0.0"
    a_expr = ["-x[0]", "0"]

    # Dirichlet BC
    uD_expr = "0"

    domain = "unit-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag


def example_3_book_stynes_2d():

    f_expr = "1.0"
    eps = 1e-3
    A_expr = [[eps,  0], [0, eps]]
    lambda_1 = eps
    lambda_expr = "10.0"
    a_expr = ["x[0]*x[0] - 1", "x[1]*x[1] - 1"]

    # Dirichlet BC
    uD_expr = "0"

    domain = "circle-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

def example_fokker_plank_paper():

    beta = 0.1
    w_plus = 2.35
    w_I = 1.9
    r = 0.3
    w_minus = 1 - r * (w_plus - 1) / (1 - r)

    f_expr = ""
    eps = beta * beta / 2
    A_expr = [[eps,  0], [0, eps]]
    W = [[w_plus - w_I,  w_plus - w_I], [w_plus - w_I, w_plus - w_I]]
    # F = nu + P(-nu + W cdot nu)
    # lambda = div F
    lambda_expr = "0.0"
    a_expr = ["x[0]*x[0] - 1", "x[1]*x[1] - 1"]

    domain = "circle-domain"
    dim = 2

    u_expr = ""
    grad_u_expr = []

    # Define whether solution is known or not
    solution_tag = "reference-solution"
    material_tag = "material-constant"

    return u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag

