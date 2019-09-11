__author__ = "Svetlana Matculevich <svetlana.v.matculevich@jyu.fi>"
__date__ = ""
__copyright__ = ""
__license__ = ""

from math import sin, cos
from dolfin import *
from numpy import *

DOLFIN_PI = 3.14159265359
def coeff_sum_term_1(alpha_n): return sin(alpha_n) / (alpha_n - sin(alpha_n) * cos(alpha_n))

def simple_example_2_2d():
    #
    p_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1])'
    p_xx_expr = '(-2) * x[1] * (1 - x[1])'
    p_yy_expr = 'x[0] * (1 - x[0]) * (-2)'

    u0_expr = '(1 - x[0]) * x[0] * (1 - x[1]) * x[1]'
    u1_expr = '(1 - x[0]) * x[0] * (1 - x[1]) * x[1]'

    g_expr = '- (' + p_xx_expr + ' + ' + p_yy_expr + ')'
    # - 2*x*(x - 1) - 2*y*(y - 1)

    f0_expr = 'alpha*x[1]*(x[0] - 1)*(x[1] - 1) - mu*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) - 4*mu*x[1]*(x[1] - 1) - lmbda*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) + alpha*x[0]*x[1]*(x[1] - 1)'
    f1_expr = 'alpha*x[0]*(x[0] - 1)*(x[1] - 1) - mu*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) - 4*mu*x[0]*(x[0] - 1) - lmbda*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) + alpha*x[0]*x[1]*(x[0] - 1)'

    #f0_expr = '- 2*x[0]*x[0] + 2*x[0]*x[1]*x[1] - 10*x[0]*x[1] + 6*x[0] - 7*x[1]*x[1] + 11*x[1] - 2'
    #f1_expr = '- 2*x[1]*x[1] + 2*x[0]*x[0]*x[1] - 10*x[0]*x[1] + 6*x[1] - 7*x[0]*x[0] + 11*x[0] - 2'


    # size of the domain
    a = 1.0
    b = 1.0
    gdim = 2
    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / a ** 2 + 1.0 / b ** 2)

    mu = 1.0        # Lame coeff. 1
    lmbda = 1.0     # Lame coeff. 2

    if gdim == 2:   # For the plaing stress
        lmbda = 2 * lmbda * mu / (lmbda + 2 * mu)

    # Givien problem data
    problem_data = dict(u_expr=[u0_expr, u1_expr],  # Exact u
                        p_expr=p_expr,  # Exact p
                        f_expr=[f0_expr, f1_expr],  # Load function
                        t_T=10.0,  # Final time
                        g_expr=g_expr,  # Source function
                        t_N_expr=['0.0', '0.0'],
                        u0_expr=['0.0', '0.0'],  # Initial u
                        p0_expr='0.0')  # Initial p

    # Domain parameters
    domain_params = dict(domain_type="rectangular-domain",
                         l_x=a,  # x_max
                         l_y=b,  # y_max
                         l_z=0.0,  # z_max = 0
                         gdim=gdim,  # dimension
                         C_FD=C_FD)

    # Material parameters
    material_params = dict(nu=1.0,  # Poinsson's ratio
                           E=1.0,  # Young modulus
                           mu=mu,  # Lame coeff. 1
                           lmbda=lmbda,  # Lame coeff. 2
                           alpha=1.0,  # Biot's constant
                           beta=1.0,  #
                           mu_f=1.0,
                           kappa=[[1, 0], [0, 1]],  # Permiability
                           kappa_inv=[[1, 0], [0, 1]],
                           min_eig_kappa=1.0)

    return problem_data, domain_params, material_params

def simple_example_2_2d_t():
    p_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * t'
    p_xx_expr = '(-2) * x[1] * (1 - x[1]) * t'
    p_yy_expr = 'x[0] * (1 - x[0]) * (-2) * t'
    p_t_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1])'

    u0_expr = '(1 - x[0]) * x[0] * (1 - x[1]) * x[1] * t'
    u0_x_expr = '(1 - 2 * x[0]) * (1 - x[1]) * x[1] * t'
    u0_xt_expr = '(1 - 2 * x[0]) * (1 - x[1]) * x[1]'

    u1_expr = '(1 - x[0]) * x[0] * (1 - x[1]) * x[1] * t'
    u1_y_expr = '(1 - x[0]) * x[0] * (1 - 2 * x[1]) * t'
    u1_yt_expr = '(1 - x[0]) * x[0] * (1 - 2 * x[1])'

    div_u_expr = u0_x_expr + ' + ' + u1_y_expr
    div_u_t_expr = u0_xt_expr + ' + ' + u1_yt_expr

    g_expr = 'beta * ' + p_t_expr + '+ alpha * ' + div_u_t_expr + '- (' + p_xx_expr + ' + ' + p_yy_expr + ')'
    # div_u_t_expr = 2*x^2*y - x^2 + 2*x*y^2 - 4*x*y + x - y^2 + y
    # g            = alpha*(2*x^2*y - x^2 + 2*x*y^2 - 4*x*y + x - y^2 + y) - 2*t*x*(x - 1) - 2*t*y*(y - 1) + beta*x*y*(x - 1)*(y - 1)
    # g            = x[0] + x[1] - 4*x[0]*x[1] + 2*x[0]*x[1]*x[1] + 2*x[0]*x[0]*x[1] - x[0]*x[0] - x[1]*x[1] - 2*t*x[0]*(x[0] - 1) - 2*t*x[1]*(x[1] - 1) + x[0]*x[1]*(x[0] - 1)*(x[1] - 1)
    # x + y - 4*x*y + 2*x*y^2 + 2*x^2*y - x^2 - y^2 - 2*t*x*(x - 1) - 2*t*y*(y - 1) + x*y*(x - 1)*(y - 1)

    f0_expr = 'alpha*t*x[1]*(x[0] - 1)*(x[1] - 1) - mu*t*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) - 4*mu*t*x[1]*(x[1] - 1) - lmbda*t*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) + alpha*t*x[0]*x[1]*(x[1] - 1)'
    f1_expr = 'alpha*t*x[0]*(x[0] - 1)*(x[1] - 1) - mu*t*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) - 4*mu*t*x[0]*(x[0] - 1) - lmbda*t*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) + alpha*t*x[0]*x[1]*(x[0] - 1)'
    # f =
    # -t*(2*x^2 - 2*x*y^2 + 10*x*y - 6*x + 7*y^2 - 11*y + 2)
    # -t*(- 2*x^2*y + 7*x^2 + 10*x*y - 11*x + 2*y^2 - 6*y + 2)

    # alph*t*x[1]*(x[0] - 1)*(x[1] - 1) - mu*t*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) - 4*mu*t*x[1]*(x[1] - 1) - lmbda*t*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) + alph*t*x[0]*x[1]*(x[1] - 1)
    # alph*t*x[0]*(x[0] - 1)*(x[1] - 1) - mu*t*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) - 4*mu*t*x[0]*(x[0] - 1) - lmbda*t*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) + alph*t*x[0]*x[1]*(x[0] - 1)


    # size of the domain
    gdim = 2
    a = 1.0
    b = 1.0
    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / a ** 2 + 1.0 / b ** 2)

    mu = 1.0        # Lame coeff. 1
    lmbda = 1.0     # Lame coeff. 2

    if gdim == 2:  # For the plane stress
        lmbda = 2 * lmbda * mu / (lmbda + 2 * mu)

    # Givien problem data
    problem_data = dict(u_expr=[u0_expr, u1_expr],  # Exact u
                        p_expr=p_expr,  # Exact p
                        f_expr=[f0_expr, f1_expr],  # Load function
                        t_T=10.0,  # Final time
                        g_expr=g_expr,  # Source function
                        t_N_expr=['0.0', '0.0'],
                        u0_expr=['0.0', '0.0'],  # Initial u
                        p0_expr='0.0')  # Initial p

    # Domain parameters
    domain_params = dict(domain_type="rectangular-domain",
                         l_x=a,  # x_max
                         l_y=b,  # y_max
                         l_z=0.0,  # z_max = 0
                         gdim=gdim,  # dimension
                         C_FD=C_FD)

    # Material parameters
    material_params = dict(nu=1.0,  # Poinsson's ratio
                           E=1.0,  # Young modulus
                           mu=mu,  # Lame coeff. 1
                           lmbda=lmbda,  # Lame coeff. 2
                           alpha=1.0,  # Biot's constant
                           beta=1.0,  #
                           mu_f=1.0,
                           kappa=[[1, 0], [0, 1]],  # Permiability
                           kappa_inv=[[1, 0], [0, 1]],
                           min_eig_kappa=1.0)

    return problem_data, domain_params, material_params

def example_1_bothetall_paper_2d_t():
    p_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * t'
    p_xx_expr = '(-2) * x[1] * (1 - x[1]) * t'
    p_yy_expr = 'x[0] * (1 - x[0]) * (-2) * t'
    p_t_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1])'

    u0_expr = '(1 - x[0]) * x[0] * (1 - x[1]) * x[1] * t'
    u0_x_expr = '(1 - 2 * x[0]) * (1 - x[1]) * x[1] * t'
    u0_xt_expr = '(1 - 2 * x[0]) * (1 - x[1]) * x[1]'

    u1_expr = '(1 - x[0]) * x[0] * (1 - x[1]) * x[1] * t'
    u1_y_expr = '(1 - x[0]) * x[0] * (1 - 2 * x[1]) * t'
    u1_yt_expr = '(1 - x[0]) * x[0] * (1 - 2 * x[1])'

    div_u_t_expr = u0_xt_expr + ' + ' + u1_yt_expr

    g_expr = 'beta * ' + p_t_expr + '+ alpha * ' + div_u_t_expr + '- (' + p_xx_expr + ' + ' + p_yy_expr + ')'
    # div_u_t_expr = 2*x^2*y - x^2 + 2*x*y^2 - 4*x*y + x - y^2 + y
    # g            = alpha*(2*x^2*y - x^2 + 2*x*y^2 - 4*x*y + x - y^2 + y) - 2*t*x*(x - 1) - 2*t*y*(y - 1) + beta*x*y*(x - 1)*(y - 1)
    # g            = x[0] + x[1] - 4*x[0]*x[1] + 2*x[0]*x[1]*x[1] + 2*x[0]*x[0]*x[1] - x[0]*x[0] - x[1]*x[1] - 2*t*x[0]*(x[0] - 1) - 2*t*x[1]*(x[1] - 1) + x[0]*x[1]*(x[0] - 1)*(x[1] - 1)
    # x + y - 4*x*y + 2*x*y^2 + 2*x^2*y - x^2 - y^2 - 2*t*x*(x - 1) - 2*t*y*(y - 1) + x*y*(x - 1)*(y - 1)

    f0_expr = 'alpha*t*x[1]*(x[0] - 1)*(x[1] - 1) - mu*t*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) - 4*mu*t*x[1]*(x[1] - 1) - lmbda*t*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) + alpha*t*x[0]*x[1]*(x[1] - 1)'
    f1_expr = 'alpha*t*x[0]*(x[0] - 1)*(x[1] - 1) - mu*t*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) - 4*mu*t*x[0]*(x[0] - 1) - lmbda*t*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) + alpha*t*x[0]*x[1]*(x[0] - 1)'
    # f =
    # -t*(2*x^2 - 2*x*y^2 + 10*x*y - 6*x + 7*y^2 - 11*y + 2)
    # -t*(- 2*x^2*y + 7*x^2 + 10*x*y - 11*x + 2*y^2 - 6*y + 2)

    # alph*t*x[1]*(x[0] - 1)*(x[1] - 1) - mu*t*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) - 4*mu*t*x[1]*(x[1] - 1) - lmbda*t*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) + alph*t*x[0]*x[1]*(x[1] - 1)
    # alph*t*x[0]*(x[0] - 1)*(x[1] - 1) - mu*t*(4*x[0]*x[1] - 4*x[1] - 2*x[0] + 2*x[1]*x[1] + 1) - 4*mu*t*x[0]*(x[0] - 1) - lmbda*t*(4*x[0]*x[1] - 2*x[1] - 4*x[0] + 2*x[0]*x[0] + 1) + alph*t*x[0]*x[1]*(x[0] - 1)

    # size of the domain
    gdim = 2
    a = 1.0
    b = 1.0
    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / a ** 2 + 1.0 / b ** 2)

    mu = 1.0        # Lame coeff. 1
    lmbda = 1.0     # Lame coeff. 2

    if gdim == 2:  # For the plaing stress
        lmbda = 2 * lmbda * mu / (lmbda + 2 * mu)

    # Givien problem data
    problem_data = dict(u_expr=[u0_expr, u1_expr],  # Exact u
                        p_expr=p_expr,  # Exact p
                        f_expr=[f0_expr, f1_expr],  # Load function
                        t_T=10.0,  # Final time
                        g_expr=g_expr,  # Source function
                        t_N_expr=['0.0', '0.0'],
                        u0_expr=['0.0', '0.0'],  # Initial u
                        p0_expr='0.0')  # Initial p

    # Domain parameters
    domain_params = dict(domain_type="rectangular-domain",
                         l_x=a,  # x_max
                         l_y=b,  # y_max
                         l_z=0.0,  # z_max = 0
                         gdim=gdim,  # dimension
                         C_FD=C_FD)

    # Material parameters
    material_params = dict(nu=1.0,  # Poinsson's ratio
                           E=1.0,  # Young modulus
                           mu=mu,  # Lame coeff. 1
                           lmbda=lmbda,  # Lame coeff. 2
                           alpha=1.0,  # Biot's constant
                           beta=1.0,  #
                           mu_f=1.0,
                           kappa=[[1, 0], [0, 1]],  # Permiability
                           kappa_inv=[[1, 0], [0, 1]],
                           min_eig_kappa=1.0)

    return problem_data, domain_params, material_params

def simple_example_2d_t():
    u0_expr = '(x[0] * x[0] + x[1] * x[1]) * t'
    u0_x_expr = '2*x[0] * t'
    u0_xt_expr = '2*x[0]'

    u1_expr = '(x[0] + x[1]) * t'
    u1_y_expr = 't'
    u1_yt_expr = '1'

    div_u_expr = u0_x_expr + ' + ' + u1_y_expr      # t*(2*x + 1)
    div_u_t_expr = u0_xt_expr + ' + ' + u1_yt_expr  # 2*x + 1

    p_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * t'
    p_x_expr = '(1 - 2 * x[0]) * x[1] * (1 - x[1]) * t'
    p_y_expr = '(1 - 2 * x[1]) * x[0] * (1 - x[0]) * t'
    p_xx_expr = '(-2) * x[1] * (1 - x[1]) * t'
    p_yy_expr = '(-2) * x[0] * (1 - x[0]) * t'
    p_t_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1])'

    g_expr = 'beta * ' + p_t_expr + '+ alpha * ' + div_u_t_expr + '- (' + p_xx_expr + ' + ' + p_yy_expr + ')'
             # bet*x*y*(x - 1)*(y - 1) + alph*(2*x + 1) - 2*t*x*(x - 1) - 2*t*y*(y - 1)

    # alph*t*y*(x - 1)*(y - 1) - 6*mu*t - 2*lmbda*t + alph*t*x*y*(y - 1)
    f_0_expr = '- mu * 6.0 * t ' \
               '- 2.0 * lmbda * t ' \
               '+ alpha * x[1] * (1 - x[1]) * (1 - 2.0 * x[0]) * t'
    # alph*t*x*(2*y - 1)*(x - 1)
    f_1_expr = 'alpha * x[0] * (1 - x[0]) * (1 - 2.0 * x[1]) * t'

    gdim = 2
    nu = 0.2           # Poinsson's ratio \in (0, 0.5)
    E = 0.594          # Young modulus
    mu = E / (2.0 * (1.0 + nu))  # Lame coeff. 1
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))  # Lame coeff. 2

    if gdim == 2:  # For the plain stress
        lmbda = 2 * lmbda * mu / (lmbda + 2 * mu)

    alpha = 1.0
    c_0 = 1e-1
    K = lmbda + 2 / 3 * mu  # skeleton bulk modulus
    K_u = K + alpha ** 2 / c_0

    min_eig_kappa = 1.0
    c_f = 1 / c_0 * min_eig_kappa * (K + 4 / 3 * mu) / (K_u + 4 / 3 * mu)
    phi_0 = 0.2
    M = 1.65e10
    beta = 1 / M + c_f * phi_0

    a = 1.0
    b = 1.0
    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / a ** 2 + 1.0 / b ** 2)

    # Givien problem data
    problem_data = dict(u_expr=[u0_expr, u1_expr],  # Exact solution
                        p_expr=p_expr,
                        f_expr=[f_0_expr, f_1_expr],  # Load function
                        t_T=10.0,  # Final time
                        g_expr=g_expr,  # Source function
                        t_N_expr=['0.0', '0.0'],
                        u0_expr=['0.0', '0.0'],
                        p0_expr='0.0')

    # Domain parameters
    domain_params = dict(domain_type="rectangular-domain",
                         l_x=a,  # a
                         l_y=b,  # b
                         l_z=0.0,
                         gdim=gdim,
                         C_FD=C_FD)

    # Material parameters
    material_params = dict(mu=mu,
                           E=E,
                           lmbda=lmbda,
                           nu=nu,  # Poinsson's ratio
                           alpha=alpha,  # Biot's constant
                           beta=beta,  #
                           kappa=[[1, 0], [0, 1]],  # permiability
                           mu_f=1.0,
                           c=c_f,
                           min_eig_kappa=1.0)

    return problem_data, domain_params, material_params


def simple_example_3d_t():
    u0_expr = '(x[0] * x[0] + x[1] * x[1] + x[2]*x[2]) * t'
    u0_x_expr = '2*x[0] * t'
    u0_xt_expr = '2*x[0]'

    u1_expr = '(x[0] + x[1] + x[2]) * t'
    u1_y_expr = 't'
    u1_yt_expr = '1'

    u2_expr = 'sin(t)'
    u2_y_expr = '0'
    u3_yt_expr = '0'

    div_u_expr = u0_x_expr + ' + ' + u1_y_expr      # t*(2*x + 1)
    div_u_t_expr = u0_xt_expr + ' + ' + u1_yt_expr  # 2*x + 1

    p_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * t'
    p_x_expr = '(1 - 2 * x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * t'
    p_y_expr = '(1 - 2 * x[1]) * x[0] * (1 - x[0]) * x[2] * (1 - x[2]) * t'
    p_z_expr = '(1 - 2 * x[2]) * x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * t'
    p_xx_expr = '(-2) * x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * t'
    p_yy_expr = '(-2) * x[0] * (1 - x[0]) * x[2] * (1 - x[2]) * t'
    p_zz_expr = '(-2) * x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * t'
    p_t_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * x[2] * (1 - x[2])'

    g_expr = 'alph*(2*x[0] + 1) ' \
             '+ 2*t*x[0]*x[1]*(x[0] - 1)*(x[1] - 1) ' \
             '+ 2*t*x[0]*x[2]*(x[0] - 1)*(x[2] - 1) ' \
             '+ 2*t*x[1]*x[2]*(x[1] - 1)*(x[2] - 1) ' \
             '- bet*x[0]*x[1]*x[2]*(x[0] - 1)*(x[1] - 1)*(x[2] - 1)'
    f_0_expr = '-2*lmbda*t - 6*mu*t - alph*t*x[0]*x[1]*x[2]*(x[1] - 1)*(x[2] - 1) ' \
               '-alph*t*x[1]*x[2]*(x[0] - 1)*(x[1] - 1)*(x[2] - 1)'
    f_1_expr = '-alph*t*x[0]*x[2]*(2*x[1] - 1)*(x[0] - 1)*(x[2] - 1)'
    f_2_expr = '-alph*t*x[0]*x[1]*(2*x[2] - 1)*(x[0] - 1)*(x[1] - 1)'

    nu = 0.2  # Poinsson's ratio
    E = 1  # Young modulus
    mu = E / (2.0 * (1.0 + nu))  # Lame coeff. 1
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))  # Lame coeff. 2

    alpha = 1.0
    c_0 = 1e-1
    K = lmbda + 2 / 3 * mu  # skeleton bulk modulus
    K_u = K + alpha ** 2 / c_0

    min_eig_kappa = 1.0
    c_f = 1 / c_0 * min_eig_kappa * (K + 4 / 3 * mu) / (K_u + 4 / 3 * mu)
    phi_0 = 0.2
    M = 1e2
    beta = 1 / M + c_f * phi_0

    a = 1.0
    b = 1.0
    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / a ** 2 + 1.0 / b ** 2)

    # Givien problem data
    problem_data = dict(u_expr=[u0_expr, u1_expr],  # Exact solution
                        p_expr=p_expr,
                        f_expr=[f_0_expr, f_1_expr, f_2_expr],  # Load function
                        t_T=10.0,  # Final time
                        g_expr=g_expr,  # Source function
                        t_N_expr=['0.0', '0.0'],
                        u0_expr=['0.0', '0.0'],
                        p0_expr='0.0')

    # Domain parameters
    domain_params = dict(domain_type="rectangular-domain",
                         l_x=a,  # a
                         l_y=b,  # b
                         l_z=0.0,
                         gdim=2,
                         C_FD=C_FD)

    # Material parameters
    material_params = dict(mu=mu,
                           E=E,
                           lmbda=lmbda,
                           nu=nu,  # Poinsson's ratio
                           alpha=alpha,  # Biot's constant
                           beta=beta,  #
                           kappa=[[1, 0], [0, 1]],  # permiability
                           kappa_inv=[[1, 0], [0, 1]],
                           mu_f=1.0,
                           c=0.465,
                           min_eig_kappa=1.0)

    return problem_data, domain_params, material_params


def mandels_problem_2d_t():
    F = 2000.0  # force
    nu_u = 0.44  # undrained poinsson's ratio
    a = 5.0  # dimentioned in x
    b = 1.0  # dimentioned in y
    c_f = 3.03 * 1e-10  # fluid comressibility (diffusivity)
    E = 1 * 1e4  # Young modulus
    nu = 0.2  # Poinsson's ratio
    c_0 = 1e-1
    alpha = 1.0
    mu = 2.475 * 1e9  # Lame coeff.
    lmbda = 1.65 * 1e9  # Lame coeff.
    K = lmbda + 2 / 3 * mu  # skeleton bulk modulus
    K_u = K + alpha ** 2 / c_0
    c_0 = 1e-1
    kappa = [[100, 0], [0, 100]]
    kappa_ = 0.01
    c_f = 1 / c_0 * kappa_ * (K + 4 / 3 * mu) / (K_u + 4 / 3 * mu)
    B = alpha / c_0 / K_u  # skemppton coeff

    # roots of (1 - nu)/(nu_u - nu) * alpha_n = tan(alpha_n)
    alpha_n = [1.35252233865, 4.6479335767, 7.8156157776952, 10.9682293796953, 14.1159175364203]

    #         2.0 *  F *  B * (1 + nu_u) / (3 * a) + sin(alpha_n) / (alpha_n - sin(alpha_n) * cos(alpha_n)) * (cos(alpha_n *x[0]/a) - cos(alpha_n)) * exp(-alpha^2 * c_f * t / a^2)
    # p_expr = "2.0 * %f * %f * (1 + %f) / (3 * %f) * " \
    #         "%f * (cos(%f *x[0]/%f) - cos(%f)) * exp( - (%f)**2 * %f * t / (%f)**2)" \
    #         %(F, B, nu_u, a, coeff_1, alpha_1, alpha_1, a, alpha_1, c_f, a)
    # Check where us the x[1] ?
    # Pressure is only changing along x[0]
    p_a_1 = '(cos(alpha_1 *x[0]/a) - cos(alpha_1))'
    p_a_2 = '(cos(alpha_2 *x[0]/a) - cos(alpha_2))'
    p_a_3 = '(cos(alpha_3 *x[0]/a) - cos(alpha_3))'
    p_a_4 = '(cos(alpha_4 *x[0]/a) - cos(alpha_4))'
    p_a_5 = '(cos(alpha_5 *x[0]/a) - cos(alpha_5))'

    p_b_1 = 'sin(alpha_1) / (alpha_1 - sin(alpha_1) * cos(alpha_1)) * exp(-alpha_1 * alpha_1 * c_f * t / a / a)'
    p_b_2 = 'sin(alpha_2) / (alpha_2 - sin(alpha_2) * cos(alpha_2)) * exp(-alpha_2 * alpha_2 * c_f * t / a / a)'
    p_b_3 = 'sin(alpha_3) / (alpha_3 - sin(alpha_3) * cos(alpha_3)) * exp(-alpha_3 * alpha_3 * c_f * t / a / a)'
    p_b_4 = 'sin(alpha_4) / (alpha_4 - sin(alpha_4) * cos(alpha_4)) * exp(-alpha_4 * alpha_4 * c_f * t / a / a)'
    p_b_5 = 'sin(alpha_5) / (alpha_5 - sin(alpha_5) * cos(alpha_5)) * exp(-alpha_5 * alpha_5 * c_f * t / a / a)'

    p_expr = '2.0 *  F *  B * (1 + nu_u) / 3.0 / a * ' \
             '(' + p_b_1 + ' * ' + p_a_1 + ' + ' \
             + p_b_2 + ' * ' + p_a_2 + ' + ' \
             + p_b_3 + ' * ' + p_a_3 + ' + ' \
             + p_b_4 + ' * ' + p_a_4 + ' + ' \
             + p_b_5 + ' * ' + p_a_5 + '' \
             + ')'

    u_a_1 = 'sin(alpha_1 *x[0]/a)'
    u_a_2 = 'sin(alpha_2 *x[0]/a)'
    u_a_3 = 'sin(alpha_3 *x[0]/a)'
    u_a_4 = 'sin(alpha_4 *x[0]/a)'
    u_a_5 = 'sin(alpha_5 *x[0]/a)'

    u_b_1 = 'sin(alpha_1)*cos(alpha_1) / (alpha_1 - sin(alpha_1) * cos(alpha_1)) * exp(-alpha_1 * alpha_1 * c_f * t / a / a)'
    u_b_2 = 'sin(alpha_2)*cos(alpha_2) / (alpha_2 - sin(alpha_2) * cos(alpha_2)) * exp(-alpha_2 * alpha_2 * c_f * t / a / a)'
    u_b_3 = 'sin(alpha_3)*cos(alpha_3) / (alpha_3 - sin(alpha_3) * cos(alpha_3)) * exp(-alpha_3 * alpha_3 * c_f * t / a / a)'
    u_b_4 = 'sin(alpha_4)*cos(alpha_4) / (alpha_4 - sin(alpha_4) * cos(alpha_4)) * exp(-alpha_4 * alpha_4 * c_f * t / a / a)'
    u_b_5 = 'sin(alpha_5)*cos(alpha_5) / (alpha_5 - sin(alpha_5) * cos(alpha_5)) * exp(-alpha_5 * alpha_5 * c_f * t / a / a)'

    u_c_1 = 'cos(alpha_1) / (alpha_1 - sin(alpha_1) * cos(alpha_1)) * exp(-alpha_1 * alpha_1 * c_f * t / a / a)'
    u_c_2 = 'cos(alpha_2) / (alpha_2 - sin(alpha_2) * cos(alpha_2)) * exp(-alpha_2 * alpha_2 * c_f * t / a / a)'
    u_c_3 = 'cos(alpha_3) / (alpha_3 - sin(alpha_3) * cos(alpha_3)) * exp(-alpha_3 * alpha_3 * c_f * t / a / a)'
    u_c_4 = 'cos(alpha_4) / (alpha_4 - sin(alpha_4) * cos(alpha_4)) * exp(-alpha_4 * alpha_4 * c_f * t / a / a)'
    u_c_5 = 'cos(alpha_5) / (alpha_5 - sin(alpha_5) * cos(alpha_5)) * exp(-alpha_5 * alpha_5 * c_f * t / a / a)'

    u_expr_0 = 'x[0] * F / mu / a * ( nu / 2 - nu_u * (' + u_b_1 + ' + ' + u_b_2 + ' + ' + u_b_3 + ' + ' + u_b_4 + ' + ' + u_b_5 + '))' \
               '+ F / mu * (' + u_c_1 + ' * ' + u_a_1 + ' + ' \
               + u_c_2 + ' * ' + u_a_2 + ' + ' \
               + u_c_3 + ' * ' + u_a_3 + ' + ' \
               + u_c_4 + ' * ' + u_a_4 + ' + ' \
               + u_c_5 + ' * ' + u_a_5 + '' \
               + ')'
    u_expr_1 = 'x[1] * F / mu / a * ( (1 - nu) / 2 - (1 - nu_u) * (' + u_b_1 + ' + ' + u_b_2 + ' + ' + u_b_3 + ' + ' + u_b_4 + ' + ' + u_b_5 + '))'

    # Givien problem data
    problem_data = dict(u_expr=[u_expr_0, u_expr_1],  # Exact solution
                        p_expr=p_expr,
                        f_expr=['0.0', '0.0'],  # Load function
                        t_T=10.0,  # Final time
                        uD_expr=['0.0', '0.0'],  # Displacement Dirichlet BC restriction
                        t_N_expr=['0.0', '-2.0 * F'],  # Force on Neumann BC restriction
                        pD_expr='0.0',  # Pressure Dirichlet BC restriction
                        g_expr='0.0',  # Source function
                        u0_expr=['0.0', '0.0'],
                        p0_expr='0.0')

    # Domain parameters
    domain_params = dict(domain_type="rectangular-domain",
                         l_x=a,  # a
                         l_y=b,  # b
                         l_z=0.0,
                         gdim=2)

    # Material parameters
    material_params = dict(F=F,
                           nu_u=nu_u,
                           B=B,
                           mu=mu,
                           lmbda=lmbda,
                           E=E,  # Young modulus
                           nu=nu,  # Poinsson's ratio
                           alpha=alpha,  # Biot's constant
                           M=1.65 * 1e10,  # Biot's modulus
                           c_f=c_f,  # fluid comressibility
                           phi_0=0.2,  # initial porosity
                           mu_f=10.0,  # fluid viscosity
                           kappa=kappa,  # permiability
                           c=0.465,  # diffusivity coeffitient
                           alpha_n=alpha_n)

    return problem_data, domain_params, material_params

def kolesov_example_2d_t():
    # Gamma_p = Gamma_1 (top) and Gamma \ Gamma_1
    p_D_expr = ['100 * sin(t)', '0']

    g_expr = '0.0'

    f_0_expr = '0'
    f_1_expr = '0.0'

    kappa = [[1, 0], [0, 1]]
    kappa_min = 1
    alpha = 1

    nu = 0.2  # Poinsson's ratio
    E = 1  # Young modulus
    mu = E / (2.0 * (1.0 + nu))  # Lame coeff. 1
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))  # Lame coeff. 2
    c_0 = 1e-1
    K = lmbda + 2 / 3 * mu  # skeleton bulk modulus
    K_u = K + alpha ** 2 / c_0
    c_f = 1 / c_0 * kappa_min * (K + 4 / 3 * mu) / (K_u + 4 / 3 * mu)

    a = 1.0
    b = 1.0

    # Givien problem data
    problem_data = dict(p_D_expr=p_D_expr,
                        f_expr=[f_0_expr, f_1_expr],  # Load function
                        t_T=DOLFIN_PI,  # Final time
                        g_expr=g_expr,  # Source function
                        t_N_expr=['0.0', '0.0'],
                        u0_expr=['0.0', '0.0'],
                        p0_expr='0.0')

    # Domain parameters
    domain_params = dict(domain_type="rectangular-domain",
                         l_x=a,  # a
                         l_y=b,  # b
                         l_z=0.0,
                         gdim=2)

    # Material parameters
    material_params = dict(mu=mu,
                           E=E,
                           lmbda=lmbda,
                           nu=nu,  # Poinsson's ratio
                           alpha=alpha,  # Biot's constant
                           M=1e2,  # Biot's modulus
                           c_f=c_f,  # fluid comressibility
                           phi_0=0.2,  # initial porosity
                           mu_f=1.0,  # fluid viscosity
                           kappa=kappa,  # permiability
                           c=0.465)

    return problem_data, domain_params, material_params


def simple_example_3_2d_t():
    u0_expr = '(x[0] * x[0] + x[1] * x[1]) * (sin(t) + 1)'
    u0_x_expr = '2*x[0] * (sin(t) + 1)'
    u0_xt_expr = '2*x[0] * cos(t)'

    u1_expr = '(x[0] + x[1]) * (sin(t) + 1)'
    u1_y_expr = 'sin(t) + 1'
    u1_yt_expr = 'cos(t)'

    div_u_expr = u0_x_expr + ' + ' + u1_y_expr      # (sin(t) + 1)*(2*x + 1)
    div_u_t_expr = u0_xt_expr + ' + ' + u1_yt_expr  # cos(t)*(2*x + 1)

    p_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * (t^2 + t + 1)'
    p_x_expr = '(1 - 2 * x[0]) * x[1] * (1 - x[1]) * (t^2 + t + 1)'
    p_y_expr = '(1 - 2 * x[1]) * x[0] * (1 - x[0]) * (t^2 + t + 1)'
    p_xx_expr = '(-2) * x[1] * (1 - x[1]) * (t^2 + t + 1)'
    p_yy_expr = '(-2) * x[0] * (1 - x[0]) * (t^2 + t + 1)'
    p_t_expr = 'x[0] * (1 - x[0]) * x[1] * (1 - x[1]) * (2*t + 1)'

    g_expr = 'beta * ' + p_t_expr + '+ alpha * ' + div_u_t_expr + '- (' + p_xx_expr + ' + ' + p_yy_expr + ')'
             # bet*x*y*(x - 1)*(y - 1)*(2*t + 1) + alph*cos(t)*(2*x + 1) - 2*(t^2 + t + 1)*x*(x - 1) - 2*(t^2 + t + 1)*y*(y - 1)

    # alph*t*y*(x - 1)*(y - 1) - 6*mu*t - 2*lmbda*t + alph*t*x*y*(y - 1)
    f_0_expr = '- mu * 6.0 * t ' \
               '- 2.0 * lmbda * t ' \
               '+ alpha * x[1] * (1 - x[1]) * (1 - 2.0 * x[0]) * t'
    # alph*t*x*(2*y - 1)*(x - 1)
    f_1_expr = 'alpha * x[0] * (1 - x[0]) * (1 - 2.0 * x[1]) * t'

    gdim = 2
    nu = 0.2           # Poinsson's ratio \in (0, 0.5)
    E = 0.594 * 1e9    # Young modulus
    mu = E / (2.0 * (1.0 + nu))  # Lame coeff. 1
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))  # Lame coeff. 2

    if gdim == 2:  # For the plain stress
        lmbda = 2 * lmbda * mu / (lmbda + 2 * mu)

    alpha = 1.0
    c_0 = 1e-1
    K = lmbda + 2 / 3 * mu  # skeleton bulk modulus
    K_u = K + alpha ** 2 / c_0

    min_eig_kappa = 1.0
    c_f = 1 / c_0 * min_eig_kappa * (K + 4 / 3 * mu) / (K_u + 4 / 3 * mu)
    phi_0 = 0.2
    M = 1.65e10
    beta = 1 / M + c_f * phi_0

    a = 1.0
    b = 1.0
    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / a ** 2 + 1.0 / b ** 2)

    # Givien problem data
    problem_data = dict(u_expr=[u0_expr, u1_expr],  # Exact solution
                        p_expr=p_expr,
                        f_expr=[f_0_expr, f_1_expr],  # Load function
                        t_T=10.0,  # Final time
                        g_expr=g_expr,  # Source function
                        t_N_expr=['0.0', '0.0'],
                        u0_expr=['0.0', '0.0'],
                        p0_expr='0.0')

    # Domain parameters
    domain_params = dict(domain_type="rectangular-domain",
                         l_x=a,  # a
                         l_y=b,  # b
                         l_z=0.0,
                         gdim=gdim,
                         C_FD=C_FD)

    # Material parameters
    material_params = dict(mu=mu,
                           E=E,
                           lmbda=lmbda,
                           nu=nu,  # Poinsson's ratio
                           alpha=alpha,  # Biot's constant
                           beta=beta,  #
                           kappa=[[1, 0], [0, 1]],  # permiability
                           mu_f=1.0,
                           c=c_f,
                           min_eig_kappa=1.0)

    return problem_data, domain_params, material_params

