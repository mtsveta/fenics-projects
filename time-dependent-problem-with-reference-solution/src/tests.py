__author__ = 'svetlana'


#-----------------------------------------------------------------------#
# test: unit domain in R^1, solution polynomial in space and time
#-----------------------------------------------------------------------#

def polynomial_in_space_time_solution_1d_t():
    u0_expr = '(1 - x[0])*x[0]'
    uD_expr = '(1 - x[0])*x[0]*(t + 1)'
    lmbd_min = 1
    A_expr = 1
    f_expr = '(1 - x[0])*x[0] + 2*(t + 1)'

    domain = "unit-domain"
    dim = 1
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

#-----------------------------------------------------------------------#
# test: unit domain in R^1, solution polynomial in space and
# trigonometric in time
#-----------------------------------------------------------------------#

def polynomial_in_space_trigonometric_in_time_solution_1d_t():
    u0_expr = '(1 - x[0])*x[0]'
    uD_expr = '(1 - x[0])*x[0]*(sin(t) + cos(t))'
    A_expr = 1
    lmbd_min = 1
    f_expr = '(1 - x[0])*x[0]*(cos(t) - sin(t)) + 2*(sin(t) + cos(t))'

    domain = "unit-domain"
    dim = 1
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

#-----------------------------------------------------------------------#
# test: unit domain in R^1, solution polynomial in space and exponential
# in time
#-----------------------------------------------------------------------#

def polynomial_in_space_exponential_in_time_1d_t():
    uD_expr = '(1 - x[0])*x[0]*exp(t)'
    u0_expr = '(1 - x[0])*x[0]'

    A_expr = 1
    lmbd_min = 1
    f_expr = '(1 - x[0])*x[0]*exp(t) + 2*exp(t)'

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return uD_expr, u0_expr, f_expr, A_expr, lmbd_min, domain, dim, T

#-----------------------------------------------------------------------#
# test: l-shape domain  in R^2, solution with singularities
#-----------------------------------------------------------------------#

def solution_with_singularities_f_trigonometric_with_resp_t_2d_t():
    uD_expr = "0"
    u0_expr = "0"
    f_expr = "sin(t)*t"

    # *exp(-0.5*(1 - (x[0] - 0.5)*(x[0] - 0.5) / 0.04 + (x[1] - 0.7)*(x[1] - 0.7) / 0.01))

    A_expr = [[1,  0], [0, 1]]
    lmbd_min = 1

    domain = "l-shape-domain"
    dim = 2
    T = 10.0

    return uD_expr, u0_expr, f_expr, A_expr, lmbd_min, domain, dim, T

#-----------------------------------------------------------------------#
# test: l-shape domain  in R^2, solution with singularities, static
#-----------------------------------------------------------------------#

def static_solution_with_singularities_2d_t():
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
    f_expr = "0.0"

    domain = "l-shape-domain"
    dim = 2
    T = 1.0

    return u_expr, grad_u_expr, f_expr, domain, dim, T

#-----------------------------------------------------------------------#
# test: l-shape domain  in R^2, solution without singularities
#-----------------------------------------------------------------------#

def solution_without_singularities_pi_shape_2d_t():
    u0_expr = "cos(x[0])*exp(x[1])*exp(t)"
    uD_expr = "cos(x[0])*exp(x[1])*exp(t)"
    f_expr = "exp(t)*(t + 1)"

    A_expr = [[10,  0], [0, 1]]
    lmbd_min = 1

    domain = "pi-shape-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

#-----------------------------------------------------------------------#
# test: unit domain in R^3, polynomial solution
#-----------------------------------------------------------------------#

def polynomial_solution_3d_t():
    u_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - x[2])*x[2]*(t*t + t + 1)'
    grad_u_x_expr = '(1 - 2*x[0])*(1 - x[1])*x[1]*(1 - x[2])*x[2]*(t*t + t + 1)'
    grad_u_y_expr = '(1 - x[0])*x[0]*(1 - 2*x[1])*(1 - x[2])*x[2]*(t*t + t + 1)'
    grad_u_z_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - 2*x[2])*(t*t + t + 1)'
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr, grad_u_z_expr]
    f_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - x[2])*x[2]*(2*t + 1) + ' \
             '2 * ((1 - x[0])*x[0]*(1 - x[1])*x[1] + (1 - x[1])*x[1]*(1 - x[2])*x[2] + ' \
             '(1 - x[2])*x[2]*(1 - x[0])*x[0])*(t*t + t + 1)'

    domain = "unit-domain"
    dim = 3
    T = 1.0

    return u_expr, grad_u_expr, f_expr, domain, dim, T

#-----------------------------------------------------------------------#
# test: l-shape domain  in R^3, solution with singularities
#-----------------------------------------------------------------------#

def solution_with_singularities_3d_t():
    u0_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "(pow(pow(x[0], 2) + pow(x[1], 2), 1.0 / 3.0) * " \
             "sin(2.0 / 3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))) + x[2])"
    uD_expr = "x[0] == 0 && x[1] == 0 ? 0.0 : " \
             "(pow(pow(x[0], 2) + pow(x[1], 2), 1.0 / 3.0) * " \
             "sin(2.0 / 3.0 * (atan2(x[1], x[0]) + (atan2(x[1], x[0]) < 0? 2*pi : 0.))) + x[2]) * " \
             "(t*t + t + 1)"
    f_expr = "0"

    A_expr = [[1,  0, 0], [0, 1, 0], [0, 0, 1]]
    lmbd_min = 1

    domain = "l-shape-domain"
    dim = 3
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

#-----------------------------------------------------------------------#
# test: l-shape domain in R^3, solution with singularities
#-----------------------------------------------------------------------#

def solution_without_singularities_3d_t():
    u_expr = "cos(x[0])*exp(x[1])*(1 - x[2])*(t*t + t + 1)"
    grad_u_x_expr = "-sin(x[0]) * exp(x[1])*(1 - x[2])*(t*t + t + 1)"
    grad_u_y_expr = "cos(x[0]) * exp(x[1])*(1 - x[2])*(t*t + t + 1)"
    grad_u_z_expr = "cos(x[0])*exp(x[1])*(-1)*(t*t + t + 1)"
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr, grad_u_z_expr]
    f_expr = "cos(x[0])*exp(x[1])*(1 - x[2])*(2*t + 1)"

    domain = "unit-domain"
    dim = 3
    T = 1.0

    return u_expr, grad_u_expr, f_expr, domain, dim, T

#-----------------------------------------------------------------------#
# test: unit domain in R^2, polynomial solution
#-----------------------------------------------------------------------#

def polynomial_solution_2d_t():

    uD_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(t*t + t + 1)'
    u0_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]'
    f_expr = '(1 - x[0])*x[0]*(1 - x[1])*x[1]*(2*t + 1) + 2 * ((1 - x[0])*x[0] + (1 - x[1])*x[1])*(t*t + t + 1)'

    A_expr = [[1,  0], [0, 1]]
    lmbd_min = 1

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

def f_one_2d_t():

    uD_expr = '0'
    u0_expr = '0'
    #u0_expr = '0'
    f_expr = '1'

    A_expr = [[10,  0], [0, 1]]
    lmbd_min = 1

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T


def f_trigonometric_with_resp_t_2d_t():

    uD_expr = '0'
    u0_expr = '0'
    #u0_expr = '0'
    f_expr = 't * sin(t) + 1'

    A_expr = [[10,  0], [0, 1]]
    lmbd_min = 1

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

def f_trigonometric_with_resp_t_lmbd_eminus1_2d_t():

    uD_expr = '0'
    u0_expr = '0'
    #u0_expr = '0'
    f_expr = 't * sin(t) + 1'

    A_expr = [[10,  0], [0, 0.1]]
    lmbd_min = 0.1

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

def f_trigonometric_with_resp_t_pi_shaped_domain_2d_t():

    uD_expr = '0'
    u0_expr = '0'
    f_expr = 't * sin(t) * sin(pi*x[0]) + t * cos(t) * sin(pi*x[1])'

    A_expr = [[1.0,  0], [0, 1.0]]
    lmbd_min = 1.0

    domain = "pi-shape-domain"
    dim = 2
    T = 2.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

def f_with_sing_with_resp_t_pi_shaped_domain_2d_t():

    uD_expr = '0'
    u0_expr = '0'
    f_expr = 't * sin(t) * exp(-100*((x[1] - 0.5)*(x[1] - 0.5) + (x[0] - 0.8)*(x[0] - 0.8)))'

    A_expr = [[1.0,  0], [0, 1.0]]
    lmbd_min = 1.0

    domain = "pi-shape-domain"
    dim = 2
    T = 2.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

def f_trigonometric_2d_t():

    uD_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    u0_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    f_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))*(cos(t) + pi*pi * (1 + 1) * (sin(t) + 1))'

    A_expr = [[1,  0], [0, 1]]
    lmbd_min = 1

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

def f_linear_2d_t():

    uD_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    u0_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    f_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))*(1 + pi*pi * (10 + 1) * (t + 1))'

    A_expr = [[10,  0], [0, 1]]
    lmbd_min = 1

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T



def f_linear_domain_with_a_whole_2d_t():

    uD_expr = '0.0'
    u0_expr = '0.0'
    f_expr = '(t + 1)'

    A_expr = [[1,  0], [0, 1]]
    lmbd_min = 1

    domain = "square-with-hole-domain"
    dim = 2
    T = 10.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T



def f_quadratic_2d_t():

    uD_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    u0_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    f_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))*(2*t + 1 + pi*pi * (10 + 1) * (t * t + t + 1))'

    A_expr = [[10,  0], [0, 1]]
    lmbd_min = 1

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

def f_quadratic_lambda_min_eminus1_2d_t():

    uD_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    u0_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    f_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))*(2*t + 1 + pi*pi * (10 + 0.1) * (t * t + t + 1))'

    A_expr = [[10,  0], [0, 0.1]]
    lmbd_min = 0.1

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T

def f_quadratic_lambda_min_eminus2_2d_t():

    uD_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    u0_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))'
    f_expr = 'sin(pi*x[0])*cos(pi*(x[1] - 0.5))*(2*t + 1 + pi*pi * (10 + 0.01) * (t * t + t + 1))'

    A_expr = [[10,  0], [0, 0.01]]
    lmbd_min = 0.01

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T
#-----------------------------------------------------------------------#
# test: unit domain in R^2, trigonometric solution
# example 8 from the AMC paper
#-----------------------------------------------------------------------#

def trigonometric_2d_t():

    u_expr = '(sin(pi*x[0])*sin(3*pi*x[1]) + sin(pi*x[1])*sin(3*pi*x[0]))*(pow(t, 3) + sin(t) + 1)'
    grad_u_x_expr = '(pi*cos(pi*x[0])*sin(3*pi*x[1]) + 3*pi*cos(3*pi*x[0])*sin(pi*x[1]))*(pow(t, 3) + sin(t) + 1)'
    grad_u_y_expr = '(pi*cos(pi*x[1])*sin(3*pi*x[0]) + 3*pi*cos(3*pi*x[1])*sin(pi*x[0]))*(pow(t, 3) + sin(t) + 1)'
    grad_u_expr = [grad_u_x_expr, grad_u_y_expr]
    f_expr = '(sin(pi*x[0])*sin(3*pi*x[1]) + sin(pi*x[1])*sin(3*pi*x[0]))*(3*pow(t, 2) + cos(t)) + ' \
        '10*pow(pi, 2)*(sin(pi*x[0])*sin(3*pi*x[1]) + sin(3*pi*x[0])*sin(pi*x[1]))*(pow(t, 3) + sin(t) + 1)'

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u_expr, grad_u_expr, f_expr, domain, dim, T

#-----------------------------------------------------------------------#
# test: unit domain in R^2, trigonometric and drasticly changing
# in time solution
# example 14 from the paper
#-----------------------------------------------------------------------#

def trigonometric_drastic_changing_in_time_2d_t():

    u0_expr = 'sin(pi*x[0])*sin(3*pi*x[1])*(t*sin(t) + 1) ' \
             '+ sin(3*pi*x[0])*sin(pi*x[1])*(cos(t) + t*sin(3*t))'
    uD_expr = 'sin(pi*x[0])*sin(3*pi*x[1])*(t*sin(t) + 1) ' \
             '+ sin(3*pi*x[0])*sin(pi*x[1])'
    f_expr = '(sin(pi*x[0])*sin(3*pi*x[1])*(sin(t) + t*cos(t)) + sin(3*pi*x[0])*sin(pi*x[1])*(sin(3*t) - sin(t) + 3*t*cos(3*t)) ) - '\
        '(-10*pow(pi, 2)*sin(pi*x[0])*sin(3*pi*x[1])*(t*sin(t) + 1) - 10*pow(pi, 2)*sin(3*pi*x[0])*sin(pi*x[1])*(cos(t) + t*sin(3*t)) )'

    A_expr = [[1,  0], [0, 1]]
    lmbd_min = 1

    domain = "unit-domain"
    dim = 2
    T = 1.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T



def f_trigonometric_circle_domain_2d_t():

    uD_expr = '0'
    u0_expr = '0'
    #u0_expr = '0'
    f_expr = 't * sin(t)*exp(-100*((x[1] - 0.2)*(x[1] - 0.2) + (x[0] + 0.3)*(x[0] + 0.3)))'

    A_expr = [[1,  0], [0, 1]]
    lmbd_min = 1

    domain = "circle-domain"
    dim = 2
    T = 10.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T


def f_trigonometric_rect_minus_circle_domain_2d_t():

    uD_expr = '0'
    u0_expr = '0'
    #u0_expr = '0'
    f_expr = 't * sin(t)*exp(-100*((x[1] - 0.8)*(x[1] - 0.8) + (x[0] + 0.8)*(x[0] + 0.8))) ' \
             '+ t*cos(t)*exp(-100*((x[1] - 0.8)*(x[1] - 0.8) + (x[0] - 0.8)*(x[0] - 0.8)))'

    A_expr = [[10,  0], [0, 1]]
    lmbd_min = 1

    domain = "rect-minus-circle-domain"
    dim = 2
    T = 10.0

    return u0_expr, uD_expr, f_expr, A_expr, lmbd_min, domain, dim, T
