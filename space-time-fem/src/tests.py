__author__ = "Svetlana Matculevich <svetlana.v.matculevich@jyu.fi>"
__date__ = ""
__copyright__ = ""
__license__ = ""


def quadratic_polynomial_solution_2d_t():

    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - x[2])*x[2]"

    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]*(1 - x[2])*x[2]"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])*(1 - x[2])*x[2]"
    du_dt_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - 2*x[2])"
    
    grad_u_expr = [du_dx_expr, du_dy_expr, du_dt_expr]

    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]*(x[2]*(1 - x[2])*x[2]"
    d2u_dy2_expr = "(-2)*(1 - x[0])*x[0]*(x[2]*(1 - x[2])*x[2]"

    # Coefficients
    A_expr = [[1,  0, 0], [0, 1, 0], [0, 0, 0]]
    min_eig_A = 1.0
    a_expr = ["0", "0", "0"]
    lambda_expr = "0.0"
    c_H = 1
    eps = min_eig_A
    
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = c_H * u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "(" + c_H + ") * (" + du_dt_expr + ")" \
             "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"
    uD_expr = u_expr
    u0_expr = u_expr

    domain = "unit-domain"
    T = 1.0

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 2 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag


def quadratic_polynomial_solution_2d_t_2():
    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*(x[2]*x[2] + x[2] + 1)"

    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]*(x[2]*x[2] + x[2] + 1)"
    du_dy_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])*(x[2]*x[2] + x[2] + 1)"
    du_dt_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 + 2*x[2])"
    grad_u_expr = [du_dx_expr, du_dy_expr, du_dt_expr]

    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]*(x[2]*x[2] + x[2] + 1)"
    d2u_dy2_expr = "(-2)*(1 - x[0])*x[0]*(x[2]*x[2] + x[2] + 1)"

    # Coefficients
    A_expr = [[1,  0, 0], [0, 1, 0], [0, 0, 0]]
    min_eig_A = 1.0
    a_expr = ["0", "0"]
    lambda_expr = "0.0"
    c_H = 1
    eps = min_eig_A
    
    divA_gradu_expr = "(" + d2u_dx2_expr + ")" + "*" + "(" + str(A_expr[0][0]) + ")" + "+" + \
                      "(" + d2u_dy2_expr + ")" + "*" + "(" + str(A_expr[1][1]) + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "(" + str(c_H) + ") * (" + du_dt_expr + ")" \
             "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")" + "+" \
                 + "(" + a_expr[1] + ")" + "*" + "(" + grad_u_expr[1] + ")"
    uD_expr = u_expr
    u0_expr = u_expr

    domain = "unit-domain"
    T = 1.0
    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 2 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

def quadratic_polynomial_solution_1d_t():

    u_expr = "(1 - x[0])*x[0]*(1 - x[1])*x[1]"
    du_dx_expr = "(1 - 2*x[0])*(1 - x[1])*x[1]"
    du_dt_expr = "(1 - x[0])*x[0]*(1 - 2*x[1])"

    grad_u_expr = [du_dx_expr, du_dt_expr]
    d2u_dx2_expr = "(-2)*(1 - x[1])*x[1]"

    # Coeffitients
    A_expr = 1
    a_expr = "0.0"
    c_H = 1.0
    min_eig_A = 1.0
    lambda_expr = "0.0"
    eps = min_eig_A

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    #f_expr = "(" + str(c_H) + ") * (" + du_dt_expr + ")" \
    #         "-" + "(" + divA_gradu_expr + ")" + "+" \
    #             + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
    #             + "(" + a_expr + ")" + "*" + "(" + grad_u_expr[0] + ")"

    f_expr = du_dt_expr + "-" + "(" + divA_gradu_expr + ")"
    uD_expr = u_expr
    u0_expr = u_expr

    domain = "unit-domain"
    T = 1.0
    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

def quadratic_polynomial_solution_1d_t_2():

    u_expr = "(1 - x[0])*x[0]*(x[1]*x[1] + x[1] + 1)"
    du_dx_expr = "(1 - 2*x[0])*(x[1]*x[1] + x[1] + 1)"
    du_dt_expr = "(1 - x[0])*x[0]*(2*x[1] + 1)"
    d2u_dx2_expr = "(-2)*(x[1]*x[1] + x[1] + 1)"
    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1
    a_expr = "0.0"
    c_H = 1.0
    min_eig_A = 1.0
    eps = 1.0
    lambda_expr = "0.0"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = u_t - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "(" + str(c_H) + ") * (" + du_dt_expr + ")" \
             "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
   
    uD_expr = u_expr
    u0_expr = u_expr

    domain = "unit-domain"
    T = 1
    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 1 + 1
    
    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

def sinusoidal_polynomial_solution_1d_t():

    u_expr = "sin(3*pi*x[0])*x[0]*exp(x[1])"
    du_dx_expr = "(sin(3*pi*x[0]) + 3*pi*x[0]*cos(3*pi*x[0]))*exp(x[1])"
    du_dt_expr = "sin(3*pi*x[0])*x[0]*exp(x[1])"
    d2u_dx2_expr = "(3*pi*cos(3*pi*x[0]) + 3*pi*(cos(3*pi*x[0]) - 3*pi*x[0]*sin(3*pi*x[0]))*exp(x[1])"
    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    c_H = 1.0
    min_eig_A = 1.0
    eps = 1.0
    lambda_expr = "0.0"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "(" + str(c_H) + ") * (" + du_dt_expr + ")" \
             "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
   
    uD_expr = u_expr
    u0_expr = u_expr

    domain = "unit-domain"
    T = 1.0
    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag


def sinusoidal_polynomial_solution_1d_t_on_long_rectangular_domain():

    u_expr = "sin(3*pi*x[0])*x[0]*(cos(x[1])*x[1] + 1.0)"
    du_dx_expr = "(sin(3*pi*x[0]) + 3*pi*x[0]*cos(3*pi*x[0]))*(cos(x[1])*x[1] + 1.0)"
    du_dt_expr = "sin(3*pi*x[0])*x[0]*(-sin(x[1])*x[1] + cos(x[1]))"
    d2u_dx2_expr = "(3*pi*cos(3*pi*x[0]) + 3*pi*(cos(3*pi*x[0]) - 3*pi*x[0]*sin(3*pi*x[0]))*(cos(x[1])*x[1] + 1.0)"
    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = ["0.0", "0.0"]
    c_H = 1.0
    min_eig_A = 1.0
    eps = 1.0
    lambda_expr = "0.0"

    # div(A * grad u)
    divA_gradu_expr = "(" + str(A_expr[0]) + ")" + "*" + "(" + d2u_dx2_expr + ")"

    # f = - div(A * grad u) + lmbd * u + a dot grad u
    f_expr = "(" + str(c_H) + ") * (" + du_dt_expr + ")" \
             "-" + "(" + divA_gradu_expr + ")" + "+" \
                 + "(" + lambda_expr + ")" + "*" + "(" + u_expr + ")" + "+" \
                 + "(" + a_expr[0] + ")" + "*" + "(" + grad_u_expr[0] + ")"
   
    uD_expr = u_expr
    u0_expr = u_expr

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 1 + 1
    domain = "unit-domain"
    T = 1.0
    
    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

#-----------------------------------------------------------------------#
# test: l-shape domain  in R^3, solution with singularities
#-----------------------------------------------------------------------#

def solution_with_singularities_2d_t():
    c_H = 1
    u_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "(pow(pow(x[0], 2) + pow(x[1], 2), 1.0 / 3.0) * " \
             "sin(2.0 / 3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.)))) * " \
             "(sin(x[2])*x[2] + 1)"
    du_dx_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "(2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    "- x[1] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))))) * (sin(x[2])*x[2] + 1)"
    du_dy_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "(2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    "+ x[1] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))))) * (sin(x[2])*x[2] + 1)"
    du_dt_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
                    "(2.0/3.0 * pow(pow(x[0], 2) + pow(x[1], 2), -2.0/3.0) * " \
                    "(x[0] * cos(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.))) " \
                    "+ x[1] * sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))))) * (sin(x[2]) + cos(x[2])*x[2])"
    grad_u_expr = [du_dx_expr, du_dy_expr, du_dt_expr]

    f_expr = "(" + str(c_H) + ") * (" + du_dt_expr + ") + " \
             "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "(pow(pow(x[0], 2) + pow(x[1], 2), 1.0/3.0) * " \
             "sin(2.0/3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0 ? 2*pi : 0.)))) * " \
             "(sin(x[2]) + cos(x[2])*x[2])"

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    c_H = 1.0
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    uD_expr = u_expr
    u0_expr = u_expr

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 1 + 1
    domain = "l-shape-domain"
    T = 1.0
    
    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag


def example_steinbach_paper_1d_t_c_1_k_1_cH_1():

    c_H = 1.0
    #u0_expr = "x[0] >= 0 && x[0] <= 0.25 ? " + \
    #         "256 * x[0] * x[0] * (1 - 4 * x[0]) * (1 - 4 * x[0]): " + \
    #          "0"
    u_1 = "1024 * ( (192 - 1*1 * pi*pi) * (1 - cos(1*pi/4)) - 24*1*pi*sin(1*pi/4) ) / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi))"

    #u_1 = "512 / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi)) * " \
    #      "(-2*pi*pi*1*1 + 42*pi*(pi*pi*1*1 - 8)*1*sin(pi*1) " \
    #       "+ (-9*pi*pi*pi*pi*1*1*1*1 + 146*pi*pi*1*1 - 384)*cos(pi*1) + 384)"
    e_1 = "exp(- (1*pi)*(1*pi)*x[1]/ " + str(c_H) + " )"
    sin_1 = "sin(1*pi*x[0])"

    u_expr = "(" + u_1 + ") * (" + e_1 + ") * (" + sin_1 + ")"
    du_dx_expr = "(" + u_1 + ") * (" + e_1 + ") * 1*pi*cos(1*pi*x[0])"
    du_dt_expr = "(" + u_1 + ") * (-(1*pi)*(1*pi)/" + str(c_H) + " ) * (" + e_1 + ") * (" + sin_1 + ")"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    uD_expr = u_expr
    u0_expr = u_expr
    f_expr = "0.0"
    T = 2.0
    domain = "rectangular-domain-1x2"

    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

def example_steinbach_paper_1d_t_c_1_k_1_cH_1_with_u0():

    c_H = 1.0
    u0_expr = "x[0] >= 0 && x[0] <= 0.25 ? " + \
             "256 * x[0] * x[0] * pow(1 - 4 * x[0], 2): " + \
              "0"
    u_1 = "1024 * ( (192 - 1*1 * pi*pi) * (1 - cos(1*pi/4)) - 24*1*pi*sin(1*pi/4) ) / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi))"

    #u_1 = "512 / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi)) * " \
    #      "(-2*pi*pi*1*1 + 42*pi*(pi*pi*1*1 - 8)*1*sin(pi*1) " \
    #       "+ (-9*pi*pi*pi*pi*1*1*1*1 + 146*pi*pi*1*1 - 384)*cos(pi*1) + 384)"
    e_1 = "exp(- (1*pi)*(1*pi)*x[1]/ " + str(c_H) + " )"
    sin_1 = "sin(1*pi*x[0])"

    u_expr = "(" + u_1 + ") * (" + e_1 + ") * (" + sin_1 + ")"
    du_dx_expr = "(" + u_1 + ") * (" + e_1 + ") * 1*pi*cos(1*pi*x[0])"
    du_dt_expr = "(" + u_1 + ") * (-(1*pi)*(1*pi)/" + str(c_H) + " ) * (" + e_1 + ") * (" + sin_1 + ")"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    uD_expr = u_expr
    f_expr = "0.0"
    T = 2.0
    domain = "rectangular-domain-1x2"

    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag



def example_steinbach_paper_1d_t_c_1_k_1_cH_100_with_u0():

    c_H = 100.0
    u0_expr = "x[0] >= 0.0 && x[0] <= 0.25 ? " + \
              "256.0 * x[0] * x[0] * (1 - 4 * x[0]) * (1 - 4 * x[0]) : " + \
              "0.0"
    u_1 = "1024 * ( (192 - 1*1 * pi*pi) * (1 - cos(1*pi/4)) - 24*1*pi*sin(1*pi/4) ) / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi))"

    #u_1 = "512 / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi)) * " \
    #      "(-2*pi*pi*1*1 + 42*pi*(pi*pi*1*1 - 8)*1*sin(pi*1) " \
    #       "+ (-9*pi*pi*pi*pi*1*1*1*1 + 146*pi*pi*1*1 - 384)*cos(pi*1) + 384)"
    e_1 = "exp(- (1*pi)*(1*pi)*x[1]/ " + str(c_H) + " )"
    sin_1 = "sin(1*pi*x[0])"

    u_expr = "(" + u_1 + ") * (" + e_1 + ") * (" + sin_1 + ")"
    du_dx_expr = "(" + u_1 + ") * (" + e_1 + ") * 1*pi*cos(1*pi*x[0])"
    du_dt_expr = "(" + u_1 + ") * (-(1*pi)*(1*pi)/" + str(c_H) + " ) * (" + e_1 + ") * (" + sin_1 + ")"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    uD_expr = u_expr
    f_expr = "0.0"
    T = 2.0
    domain = "rectangular-domain-1x2"

    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag


'''
def example_steinbach_paper_1d_t():

    c_H = 1.0
    eps = 1.0
    u0_expr = "x[0] >= 0 && x[0] <= 0.25 ? " + \
              "256 * x[0] * x[0] * pow(1 - 4 * x[0], 2): " + \
              "0"
    u_1 = "1024 * ( (192 - 1*1 * pi*pi) * (1 - cos(1*pi/4)) - 24*1*pi*sin(1*pi/4) ) / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi))"
    u_2 = "1024 * ( (192 - 2*2 * pi*pi) * (1 - cos(2*pi/4)) - 24*2*pi*sin(2*pi/4) ) / ((2*pi)*(2*pi)*(2*pi)*(2*pi)*(2*pi))"
    u_3 = "1024 * ( (192 - 3*3 * pi*pi) * (1 - cos(3*pi/4)) - 24*3*pi*sin(3*pi/4) ) / ((3*pi)*(3*pi)*(3*pi)*(3*pi)*(3*pi))"
    u_4 = "1024 * ( (192 - 4*4 * pi*pi) * (1 - cos(4*pi/4)) - 24*4*pi*sin(4*pi/4) ) / ((4*pi)*(4*pi)*(4*pi)*(4*pi)*(4*pi))"

    e_1 = "exp(- (1*pi)*(1*pi)*x[1]/ " + str(c_H) + " )"
    e_2 = "exp(- (2*pi)*(2*pi)*x[1]/ " + str(c_H) + " )"
    e_3 = "exp(- (3*pi)*(3*pi)*x[1]/ " + str(c_H) + " )"
    e_4 = "exp(- (4*pi)*(4*pi)*x[1]/ " + str(c_H) + " )"

    sin_1 = "sin(1*pi*x[0])"
    sin_2 = "sin(2*pi*x[0])"
    sin_3 = "sin(3*pi*x[0])"
    sin_4 = "sin(4*pi*x[0])"

    u_expr = "(" + u_1 + ") * (" + e_1 + ") * (" + sin_1 + ") + " \
             "(" + u_2 + ") * (" + e_2 + ") * (" + sin_2 + ") + " \
             "(" + u_3 + ") * (" + e_3 + ") * (" + sin_3 + ") + " \
             "(" + u_4 + ") * (" + e_4 + ") * (" + sin_4 + ")"

    du_dx_expr = "(" + u_1 + ") * (" + e_1 + ") * 1*pi*cos(1*pi*x[0]) + " \
                 "(" + u_2 + ") * (" + e_2 + ") * 2*pi*cos(2*pi*x[0]) + " \
                 "(" + u_3 + ") * (" + e_3 + ") * 3*pi*cos(3*pi*x[0]) + " \
                 "(" + u_4 + ") * (" + e_4 + ") * 4*pi*cos(4*pi*x[0])"
    du_dt_expr = "(" + u_1 + ") * (-(1*pi)*(1*pi)/" + str(c_H) + " ) * (" + e_1 + ") * (" + sin_1 + ") + " \
                 "(" + u_2 + ") * (-(2*pi)*(2*pi)/" + str(c_H) + " ) * (" + e_2 + ") * (" + sin_2 + ") + " \
                 "(" + u_3 + ") * (-(3*pi)*(3*pi)/" + str(c_H) + " ) * (" + e_3 + ") * (" + sin_3 + ") + " \
                 "(" + u_4 + ") * (-(4*pi)*(4*pi)/" + str(c_H) + " ) * (" + e_4 + ") * (" + sin_4 + ")"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    f_expr = "0.0"
    T = 2.0
    domain = "rectangular-domain-1x2"

    return u_expr, grad_u_expr, f_expr, u0_expr, c_H, eps, T, domain

def example_steinbach_paper_1d_t_c_10():

    c_H = 10.0
    eps = 1.0
    u0_expr = "x[0] >= 0 && x[0] <= 0.25 ? " + \
              "256 * x[0] * x[0] * pow(1 - 4 * x[0], 2): " + \
              "0"
    u_1 = "1024 * ( (192 - 1*1 * pi*pi) * (1 - cos(1*pi/4)) - 24*1*pi*sin(1*pi/4) ) / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi))"
    u_2 = "1024 * ( (192 - 2*2 * pi*pi) * (1 - cos(2*pi/4)) - 24*2*pi*sin(2*pi/4) ) / ((2*pi)*(2*pi)*(2*pi)*(2*pi)*(2*pi))"
    u_3 = "1024 * ( (192 - 3*3 * pi*pi) * (1 - cos(3*pi/4)) - 24*3*pi*sin(3*pi/4) ) / ((3*pi)*(3*pi)*(3*pi)*(3*pi)*(3*pi))"
    u_4 = "1024 * ( (192 - 4*4 * pi*pi) * (1 - cos(4*pi/4)) - 24*4*pi*sin(4*pi/4) ) / ((4*pi)*(4*pi)*(4*pi)*(4*pi)*(4*pi))"

    e_1 = "exp(- (1*pi)*(1*pi)*x[1]/ " + str(c_H) + " )"
    e_2 = "exp(- (2*pi)*(2*pi)*x[1]/ " + str(c_H) + " )"
    e_3 = "exp(- (3*pi)*(3*pi)*x[1]/ " + str(c_H) + " )"
    e_4 = "exp(- (4*pi)*(4*pi)*x[1]/ " + str(c_H) + " )"

    sin_1 = "sin(1*pi*x[0])"
    sin_2 = "sin(2*pi*x[0])"
    sin_3 = "sin(3*pi*x[0])"
    sin_4 = "sin(4*pi*x[0])"

    u_expr = "(" + u_1 + ") * (" + e_1 + ") * (" + sin_1 + ") + " \
             "(" + u_2 + ") * (" + e_2 + ") * (" + sin_2 + ") + " \
             "(" + u_3 + ") * (" + e_3 + ") * (" + sin_3 + ") + " \
             "(" + u_4 + ") * (" + e_4 + ") * (" + sin_4 + ")"

    du_dx_expr = "(" + u_1 + ") * (" + e_1 + ") * 1*pi*cos(1*pi*x[0]) + " \
                 "(" + u_2 + ") * (" + e_2 + ") * 2*pi*cos(2*pi*x[0]) + " \
                 "(" + u_3 + ") * (" + e_3 + ") * 3*pi*cos(3*pi*x[0]) + " \
                 "(" + u_4 + ") * (" + e_4 + ") * 4*pi*cos(4*pi*x[0])"
    du_dt_expr = "(" + u_1 + ") * (-(1*pi)*(1*pi)/" + str(c_H) + " ) * (" + e_1 + ") * (" + sin_1 + ") + " \
                 "(" + u_2 + ") * (-(2*pi)*(2*pi)/" + str(c_H) + " ) * (" + e_2 + ") * (" + sin_2 + ") + " \
                 "(" + u_3 + ") * (-(3*pi)*(3*pi)/" + str(c_H) + " ) * (" + e_3 + ") * (" + sin_3 + ") + " \
                 "(" + u_4 + ") * (-(4*pi)*(4*pi)/" + str(c_H) + " ) * (" + e_4 + ") * (" + sin_4 + ")"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    f_expr = "0.0"
    T = 2.0
    domain = "rectangular-domain-1x2"

    return u_expr, grad_u_expr, f_expr, u0_expr, c_H, eps, T, domain

def example_steinbach_paper_1d_t_c_100():

    c_H = 100.0
    eps = 1.0
    u0_expr = "(x[0] >= 0.0 && x[0] <= 0.25 && x[1] >= 0.0 && x[1] <= 2.0) ? " + \
              "256.0 * x[0] * x[0] * (1 - 4 * x[0]) * (1 - 4 * x[0]) : " + \
              "0.0"
    u_1 = "1024 * ( (192 - 1*1 * pi*pi) * (1 - cos(1*pi/4)) - 24*1*pi*sin(1*pi/4) ) / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi))"
    u_2 = "1024 * ( (192 - 2*2 * pi*pi) * (1 - cos(2*pi/4)) - 24*2*pi*sin(2*pi/4) ) / ((2*pi)*(2*pi)*(2*pi)*(2*pi)*(2*pi))"
    u_3 = "1024 * ( (192 - 3*3 * pi*pi) * (1 - cos(3*pi/4)) - 24*3*pi*sin(3*pi/4) ) / ((3*pi)*(3*pi)*(3*pi)*(3*pi)*(3*pi))"
    u_4 = "1024 * ( (192 - 4*4 * pi*pi) * (1 - cos(4*pi/4)) - 24*4*pi*sin(4*pi/4) ) / ((4*pi)*(4*pi)*(4*pi)*(4*pi)*(4*pi))"
    u_5 = "1024 * ( (192 - 5*5 * pi*pi) * (1 - cos(5*pi/4)) - 24*5*pi*sin(5*pi/4) ) / ((5*pi)*(5*pi)*(5*pi)*(5*pi)*(5*pi))"
    u_6 = "1024 * ( (192 - 6*6 * pi*pi) * (1 - cos(6*pi/4)) - 24*6*pi*sin(6*pi/4) ) / ((6*pi)*(6*pi)*(6*pi)*(6*pi)*(6*pi))"

    e_1 = "exp(- (1*pi)*(1*pi)*x[1]/ " + str(c_H) + " )"
    e_2 = "exp(- (2*pi)*(2*pi)*x[1]/ " + str(c_H) + " )"
    e_3 = "exp(- (3*pi)*(3*pi)*x[1]/ " + str(c_H) + " )"
    e_4 = "exp(- (4*pi)*(4*pi)*x[1]/ " + str(c_H) + " )"
    e_5 = "exp(- (5*pi)*(5*pi)*x[1]/ " + str(c_H) + " )"
    e_6 = "exp(- (6*pi)*(6*pi)*x[1]/ " + str(c_H) + " )"

    sin_1 = "sin(1*pi*x[0])"
    sin_2 = "sin(2*pi*x[0])"
    sin_3 = "sin(3*pi*x[0])"
    sin_4 = "sin(4*pi*x[0])"
    sin_5 = "sin(5*pi*x[0])"
    sin_6 = "sin(6*pi*x[0])"

    u_expr = "(" + u_1 + ") * (" + e_1 + ") * (" + sin_1 + ") + " \
             "(" + u_2 + ") * (" + e_2 + ") * (" + sin_2 + ") + " \
             "(" + u_3 + ") * (" + e_3 + ") * (" + sin_3 + ") + " \
             "(" + u_4 + ") * (" + e_4 + ") * (" + sin_4 + ") + " \
             "(" + u_5 + ") * (" + e_5 + ") * (" + sin_5 + ") + " \
             "(" + u_6 + ") * (" + e_6 + ") * (" + sin_6 + ")"
    du_dx_expr = "(" + u_1 + ") * (" + e_1 + ") * 1*pi*cos(1*pi*x[0]) + " \
                 "(" + u_2 + ") * (" + e_2 + ") * 2*pi*cos(2*pi*x[0]) + " \
                 "(" + u_3 + ") * (" + e_3 + ") * 3*pi*cos(3*pi*x[0]) + " \
                 "(" + u_4 + ") * (" + e_4 + ") * 4*pi*cos(4*pi*x[0]) + " \
                 "(" + u_5 + ") * (" + e_5 + ") * 5*pi*cos(5*pi*x[0]) + " \
                 "(" + u_6 + ") * (" + e_6 + ") * 6*pi*cos(6*pi*x[0])"

    du_dt_expr = "(" + u_1 + ") * (-(1*pi)*(1*pi)/" + str(c_H) + " ) * (" + e_1 + ") * (" + sin_1 + ") + " \
                 "(" + u_2 + ") * (-(2*pi)*(2*pi)/" + str(c_H) + " ) * (" + e_2 + ") * (" + sin_2 + ") + " \
                 "(" + u_3 + ") * (-(3*pi)*(3*pi)/" + str(c_H) + " ) * (" + e_3 + ") * (" + sin_3 + ") + " \
                 "(" + u_4 + ") * (-(4*pi)*(4*pi)/" + str(c_H) + " ) * (" + e_4 + ") * (" + sin_4 + ") + " \
                 "(" + u_5 + ") * (-(5*pi)*(5*pi)/" + str(c_H) + " ) * (" + e_5 + ") * (" + sin_5 + ") + " \
                 "(" + u_6 + ") * (-(6*pi)*(6*pi)/" + str(c_H) + " ) * (" + e_6 + ") * (" + sin_6 + ")"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    f_expr = "0.0"
    domain = "rectangular-domain-1x2"

    T = 2.0
    return u_expr, grad_u_expr, f_expr, u0_expr, c_H, eps, T, domain

def example_steinbach_paper_1d_t_c_100_without_u0():

    c_H = 100.0
    eps = 1.0
    u_1 = "1024 * ( (192 - 1*1 * pi*pi) * (1 - cos(1*pi/4)) - 24*1*pi*sin(1*pi/4) ) / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi))"
    u_2 = "1024 * ( (192 - 2*2 * pi*pi) * (1 - cos(2*pi/4)) - 24*2*pi*sin(2*pi/4) ) / ((2*pi)*(2*pi)*(2*pi)*(2*pi)*(2*pi))"
    u_3 = "1024 * ( (192 - 3*3 * pi*pi) * (1 - cos(3*pi/4)) - 24*3*pi*sin(3*pi/4) ) / ((3*pi)*(3*pi)*(3*pi)*(3*pi)*(3*pi))"
    u_4 = "1024 * ( (192 - 4*4 * pi*pi) * (1 - cos(4*pi/4)) - 24*4*pi*sin(4*pi/4) ) / ((4*pi)*(4*pi)*(4*pi)*(4*pi)*(4*pi))"

    e_1 = "exp(- (1*pi)*(1*pi)*x[1]/ " + str(c_H) + " )"
    e_2 = "exp(- (2*pi)*(2*pi)*x[1]/ " + str(c_H) + " )"
    e_3 = "exp(- (3*pi)*(3*pi)*x[1]/ " + str(c_H) + " )"
    e_4 = "exp(- (4*pi)*(4*pi)*x[1]/ " + str(c_H) + " )"

    sin_1 = "sin(1*pi*x[0])"
    sin_2 = "sin(2*pi*x[0])"
    sin_3 = "sin(3*pi*x[0])"
    sin_4 = "sin(4*pi*x[0])"

    u_expr = "(" + u_1 + ") * (" + e_1 + ") * (" + sin_1 + ") + " \
             "(" + u_2 + ") * (" + e_2 + ") * (" + sin_2 + ") + " \
             "(" + u_3 + ") * (" + e_3 + ") * (" + sin_3 + ") + " \
             "(" + u_4 + ") * (" + e_4 + ") * (" + sin_4 + ")"

    du_dx_expr = "(" + u_1 + ") * (" + e_1 + ") * 1*pi*cos(1*pi*x[0]) + " \
                 "(" + u_2 + ") * (" + e_2 + ") * 2*pi*cos(2*pi*x[0]) + " \
                 "(" + u_3 + ") * (" + e_3 + ") * 3*pi*cos(3*pi*x[0]) + " \
                 "(" + u_4 + ") * (" + e_4 + ") * 4*pi*cos(4*pi*x[0])"
    du_dt_expr = "(" + u_1 + ") * (-(1*pi)*(1*pi)/" + str(c_H) + " ) * (" + e_1 + ") * (" + sin_1 + ") + " \
                 "(" + u_2 + ") * (-(2*pi)*(2*pi)/" + str(c_H) + " ) * (" + e_2 + ") * (" + sin_2 + ") + " \
                 "(" + u_3 + ") * (-(3*pi)*(3*pi)/" + str(c_H) + " ) * (" + e_3 + ") * (" + sin_3 + ") + " \
                 "(" + u_4 + ") * (-(4*pi)*(4*pi)/" + str(c_H) + " ) * (" + e_4 + ") * (" + sin_4 + ")"

    u0_expr = u_expr

    grad_u_expr = [du_dx_expr, du_dt_expr]

    f_expr = "0.0"
    domain = "rectangular-domain-1x2"

    T = 2.0
    return u_expr, grad_u_expr, f_expr, u0_expr, c_H, eps, T, domain
'''

# Paul's online notes
# Solving heat equation
def example_tutorial_example_2():
    c_H = 1.0
    u0_expr = "6*sin(pi*x[0])"
    u_expr = "6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dx_expr = "6*pi*cos(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dt_expr = "(- pi*pi/ " + str(c_H) + " )*6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    f_expr = "0.0"
    domain = "rectangular-domain-1x2"
    uD_expr = u_expr

    T = 2.0
    dim = 1 + 1

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"


    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag


def example_tutorial_example_2_cH_1():

    c_H = 1.0

    u0_expr = "6*sin(pi*x[0])"
    u_expr = "6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dx_expr = "6*pi*cos(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dt_expr = "(- pi*pi/ " + str(c_H) + " )*6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    f_expr = "0.0"
    domain = "rectangular-domain-1x2"
    uD_expr = u_expr

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"

    T = 2.0
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

def example_tutorial_example_2_cH_10():

    c_H = 10.0

    u0_expr = "6*sin(pi*x[0])"
    u_expr = "6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dx_expr = "6*pi*cos(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dt_expr = "(- pi*pi/ " + str(c_H) + " )*6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    f_expr = "0.0"
    domain = "rectangular-domain-1x2"
    uD_expr = u_expr

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"

    T = 2.0
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

def example_tutorial_example_2_cH_20():

    c_H = 20.0
    u0_expr = "6*sin(pi*x[0])"
    u_expr = "6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dx_expr = "6*pi*cos(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dt_expr = "(- pi*pi/ " + str(c_H) + " )*6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    f_expr = "0.0"
    domain = "rectangular-domain-1x2"
    uD_expr = u_expr

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"

    T = 2.0
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

def example_tutorial_example_2_cH_100():

    c_H = 100.0
    u0_expr = "6*sin(pi*x[0])"
    u_expr = "6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dx_expr = "6*pi*cos(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"
    du_dt_expr = "(- pi*pi/ " + str(c_H) + " )*6*sin(pi*x[0])*exp(- pi*pi*x[1]/ " + str(c_H) + " )"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    f_expr = "0.0"
    domain = "rectangular-domain-1x2"
    uD_expr = u_expr

    # Define whether solution is known or not
    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"

    T = 2.0
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

def example_steinbach_paper_1d_t_c_1_k_1_cH_10():

    c_H = 10.0
    #u0_expr = "x[0] >= 0 && x[0] <= 0.25 ? " + \
    #         "256 * x[0] * x[0] * (1 - 4 * x[0]) * (1 - 4 * x[0]): " + \
    #          "0"
    u_1 = "1024 * ( (192 - 1*1 * pi*pi) * (1 - cos(1*pi/4)) - 24*1*pi*sin(1*pi/4) ) / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi))"

    #u_1 = "512 / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi)) * " \
    #      "(-2*pi*pi*1*1 + 42*pi*(pi*pi*1*1 - 8)*1*sin(pi*1) " \
    #       "+ (-9*pi*pi*pi*pi*1*1*1*1 + 146*pi*pi*1*1 - 384)*cos(pi*1) + 384)"
    e_1 = "exp(- (1*pi)*(1*pi)*x[1]/ " + str(c_H) + " )"
    sin_1 = "sin(1*pi*x[0])"

    u_expr = "(" + u_1 + ") * (" + e_1 + ") * (" + sin_1 + ")"
    du_dx_expr = "(" + u_1 + ") * (" + e_1 + ") * 1*pi*cos(1*pi*x[0])"
    du_dt_expr = "(" + u_1 + ") * (-(1*pi)*(1*pi)/" + str(c_H) + " ) * (" + e_1 + ") * (" + sin_1 + ")"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    # Coeffitients
    A_expr = 1.0
    a_expr = "0.0"
    min_eig_A = 1.0
    eps = min_eig_A
    lambda_expr = "0.0"

    uD_expr = u_expr
    u0_expr = u_expr
    f_expr = "0.0"
    T = 2.0
    domain = "rectangular-domain-1x2"

    solution_tag = "predefined-solution"
    material_tag = "material-constant"
    pde_tag = "heat-pde"
    dim = 1 + 1

    return u_expr, grad_u_expr, f_expr, \
           A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
           u0_expr, uD_expr, T, domain, dim, \
           solution_tag, material_tag, pde_tag

def example_steinbach_paper_1d_t_c_100():

    c_H = 100.0
    eps = 1.0
    u0_expr = "x[0] >= 0.0 && x[0] <= 0.25 ? " + \
              "256.0 * x[0] * x[0] * (1 - 4 * x[0]) * (1 - 4 * x[0]) : " + \
              "0.0"
    u_1 = "1024 * ( (192 - 1*1 * pi*pi) * (1 - cos(1*pi/4)) - 24*1*pi*sin(1*pi/4) ) / ((1*pi)*(1*pi)*(1*pi)*(1*pi)*(1*pi))"
    u_2 = "1024 * ( (192 - 2*2 * pi*pi) * (1 - cos(2*pi/4)) - 24*2*pi*sin(2*pi/4) ) / ((2*pi)*(2*pi)*(2*pi)*(2*pi)*(2*pi))"
    u_3 = "1024 * ( (192 - 3*3 * pi*pi) * (1 - cos(3*pi/4)) - 24*3*pi*sin(3*pi/4) ) / ((3*pi)*(3*pi)*(3*pi)*(3*pi)*(3*pi))"
    u_4 = "1024 * ( (192 - 4*4 * pi*pi) * (1 - cos(4*pi/4)) - 24*4*pi*sin(4*pi/4) ) / ((4*pi)*(4*pi)*(4*pi)*(4*pi)*(4*pi))"
    u_5 = "1024 * ( (192 - 5*5 * pi*pi) * (1 - cos(5*pi/4)) - 24*5*pi*sin(5*pi/4) ) / ((5*pi)*(5*pi)*(5*pi)*(5*pi)*(5*pi))"
    u_6 = "1024 * ( (192 - 6*6 * pi*pi) * (1 - cos(6*pi/4)) - 24*6*pi*sin(6*pi/4) ) / ((6*pi)*(6*pi)*(6*pi)*(6*pi)*(6*pi))"

    e_1 = "exp(- (1*pi)*(1*pi)*x[1]/ " + str(c_H) + " )"
    e_2 = "exp(- (2*pi)*(2*pi)*x[1]/ " + str(c_H) + " )"
    e_3 = "exp(- (3*pi)*(3*pi)*x[1]/ " + str(c_H) + " )"
    e_4 = "exp(- (4*pi)*(4*pi)*x[1]/ " + str(c_H) + " )"
    e_5 = "exp(- (5*pi)*(5*pi)*x[1]/ " + str(c_H) + " )"
    e_6 = "exp(- (6*pi)*(6*pi)*x[1]/ " + str(c_H) + " )"

    sin_1 = "sin(1*pi*x[0])"
    sin_2 = "sin(2*pi*x[0])"
    sin_3 = "sin(3*pi*x[0])"
    sin_4 = "sin(4*pi*x[0])"
    sin_5 = "sin(5*pi*x[0])"
    sin_6 = "sin(6*pi*x[0])"

    u_expr = "(" + u_1 + ") * (" + e_1 + ") * (" + sin_1 + ") + " \
             "(" + u_2 + ") * (" + e_2 + ") * (" + sin_2 + ") + " \
             "(" + u_3 + ") * (" + e_3 + ") * (" + sin_3 + ") + " \
             "(" + u_4 + ") * (" + e_4 + ") * (" + sin_4 + ") + " \
             "(" + u_5 + ") * (" + e_5 + ") * (" + sin_5 + ") + " \
             "(" + u_6 + ") * (" + e_6 + ") * (" + sin_6 + ")"
    du_dx_expr = "(" + u_1 + ") * (" + e_1 + ") * 1*pi*cos(1*pi*x[0]) + " \
                 "(" + u_2 + ") * (" + e_2 + ") * 2*pi*cos(2*pi*x[0]) + " \
                 "(" + u_3 + ") * (" + e_3 + ") * 3*pi*cos(3*pi*x[0]) + " \
                 "(" + u_4 + ") * (" + e_4 + ") * 4*pi*cos(4*pi*x[0]) + " \
                 "(" + u_5 + ") * (" + e_5 + ") * 5*pi*cos(5*pi*x[0]) + " \
                 "(" + u_6 + ") * (" + e_6 + ") * 6*pi*cos(6*pi*x[0])"

    du_dt_expr = "(" + u_1 + ") * (-(1*pi)*(1*pi)/" + str(c_H) + " ) * (" + e_1 + ") * (" + sin_1 + ") + " \
                 "(" + u_2 + ") * (-(2*pi)*(2*pi)/" + str(c_H) + " ) * (" + e_2 + ") * (" + sin_2 + ") + " \
                 "(" + u_3 + ") * (-(3*pi)*(3*pi)/" + str(c_H) + " ) * (" + e_3 + ") * (" + sin_3 + ") + " \
                 "(" + u_4 + ") * (-(4*pi)*(4*pi)/" + str(c_H) + " ) * (" + e_4 + ") * (" + sin_4 + ") + " \
                 "(" + u_5 + ") * (-(5*pi)*(5*pi)/" + str(c_H) + " ) * (" + e_5 + ") * (" + sin_5 + ") + " \
                 "(" + u_6 + ") * (-(6*pi)*(6*pi)/" + str(c_H) + " ) * (" + e_6 + ") * (" + sin_6 + ")"

    grad_u_expr = [du_dx_expr, du_dt_expr]

    f_expr = "0.0"
    domain = "rectangular-domain-1x2"

    T = 2.0
    return u_expr, grad_u_expr, f_expr, u0_expr, c_H, eps, T, domain
