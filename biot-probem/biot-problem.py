__author__ = "Svetlana Matculevich <svetlana.v.matculevich@jyu.fi>"
__date__ = "2016-04-20"
__copyright__ = "Copyright (C) 2016 Svetlana Matculevich"
__license__ = ""

from dolfin import *
from numpy import *
from time import process_time

# Customer libraries
import tests, problem, postprocess
#from problem import epsilon, sigma, inv_sigma
#from dolfin.cpp.la import list_linear_solver_methods, list_krylov_solver_preconditioners
#from matplotlib import interactive

# https://fenicsproject.org/olddocs/dolfin/1.3.0/python/programmers-reference/fem/solving/solve.html
# https://fenicsproject.org/pub/tutorial/sphinx1/._ftut1006.html#ftut-app-solver-prec

parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['cpp_optimize'] = True

# print(parameters.linear_algebra_backend)
# list_linear_solver_methods()
# list_krylov_solver_preconditioners()

# Class to generate the errors and their estimates
class ErrorControl():

    # Init variable of ErrorControl class
    def __init__(self, mesh, VVe, Qe, gdim):
        """
        :param mesh: mesh
        :param Ve: functional space of exact displacement
        :param Qe: functional space of exact pressure
        :param gdim: dimension
        """
        self.mesh = mesh
        self.VVe = VVe
        self.Qe = Qe
        self.gdim = gdim

    def get_error_array(self, num):
        return zeros(num)

    # Update the mesh
    def update(self, mesh, VVe, Qe):
        """
        :param mesh: mesh
        :param Ve: functional space of exact displacement
        :param Qe: functional space of exact pressure
        """
        # Update the mesh
        self.mesh = mesh
        # Update the functional spaces
        self.VVe = VVe
        self.Qe = Qe

    def error_pressure(self, error_array, errorl2_array, i, p, norm_p_accum, n, funcs):
        """
        :param p: approximate pressure
        :param pe: exact pressure
        :param mu_f: fluid ...
        :param kappa: permeability tensor
        :param beta: (1/M + phi_0 * c_f) fluid characteristic
        :param Ve: functional space of exact solution
        :param mesh: mesh
        :param dim: dimension of the problem

        :return: error-norm between p and ue
        """
        # Interpolate exact and approximate solution to the functional space of exact solution
        p_ve = interpolate(p, self.Qe)
        p_exact_ve = interpolate(funcs["p_e"], self.Qe)
        e_p = p_ve - p_exact_ve

        # Define variational form of the error
        var_grad_e_p = inner(funcs["kappa"] * grad(e_p), grad(e_p))
        var_e_p      = inner(funcs["beta"] * e_p, e_p)

        # Assembling variational form of the error
        norm_grad_e_p = assemble(var_grad_e_p * dx)
        norm_e_p = assemble(var_e_p * dx)

        if test_params["error_format"] == "relative_error":
            norm_p = assemble((inner(funcs["kappa"] * grad(p_exact_ve), grad(p_exact_ve))
                               + inner(funcs["beta"] * p_exact_ve, p_exact_ve)) * dx)
        elif test_params["error_format"] == "absolute_error":
            norm_p = 1

        error_array[i-1] = norm_grad_e_p + norm_e_p
        errorl2_array[i-1] = norm_e_p

        if test_params["full_documentation"]:
            # print('%--------------------')
            # print('% Error in pressure:')
            # print('%--------------------')
            print("||| e_p |||^2    = %.4e" % (error_array[i-1] / norm_p))
            # print("|| grad e_p ||^2 = %.4e" % norm_grad_e_p)
            print("|| e_p ||^2_beta = %.4e" % (errorl2_array[i-1] / norm_p))
            # print("||| p |||^2    = %.4e" % norm_p)
            # print("|| grad e_p ||^2 = %.4e" % norm_grad_e_p)

        if n == 0:  norm_p_accum[n] = norm_p
        else:   norm_p_accum[n] = norm_p_accum[n - 1] + norm_p

        return error_array, errorl2_array, var_grad_e_p, var_e_p, norm_p, norm_p_accum

    def error_displacement(self, eu_array, eudiv_array, i, u, norm_u_accum, n, funcs, mu, lmbda):
        """
        :param eu_array: array with energy error
        :param eudiv_array: array with div error
        :param u: approximate solution
        :param ue: exact solution
        :param mu: first Lame coefficient
        :param lmbda: Lame coefficient
        :return: error-norm between u and ue
        """
        # Interpolate exact and approximate solution to the functional space of exact solution
        u_ve = interpolate(u, self.VVe)
        u_exact_ve = interpolate(funcs["u_e"], self.VVe)
        e = u_ve - u_exact_ve

        # Define variational form of the error
        epsilon = lambda v: 0.5 * (nabla_grad(v) + nabla_grad(v).T)
        sigma = lambda v: 2 * mu * epsilon(v) + lmbda * div(v) * Identity(self.gdim)

        var_energy_e_u = tr(inner(sigma(e), epsilon(e)))
        var_strain_e_u = tr(inner(epsilon(e), epsilon(e)))
        var_div_e_u    = inner(div(e), div(e)) #tr(epsilon(e))**2 #

        # Assembling variational form of the error
        energy_e_u = assemble(var_energy_e_u * dx)
        div_e_u    = assemble(var_div_e_u * dx)
        strain_e_u = assemble(var_strain_e_u * dx)

        if test_params["error_format"] == "relative_error":
            norm_u = assemble(tr(inner(sigma(u_exact_ve), epsilon(u_exact_ve))) * dx)
            #norm_u = assemble((2.0 * mu * inner(epsilon(u_exact_ve), epsilon(u_exact_ve)) \
            #                   + lmbda * tr(epsilon((u_exact_ve)))**2) * dx)
        elif test_params["error_format"] == "absolute_error":
            norm_u = 1

        eu_array[i-1] = energy_e_u
        eudiv_array[i-1] = div_e_u

        if test_params["full_documentation"]:
            # print('%------------------------')
            # print('% Error in displacement:')
            # print('%------------------------')
            print("||| e_u |||^2   = %.4e" % (eu_array[i-1] / norm_u))
            # print("||  eps(e_u) ||^2 = %.4e" % strain_e_u)
            print("|| div e_p ||^2 = %.4e" % (eudiv_array[i-1] / norm_u))
            # print("||| u |||^2   = %.4e" % norm_u_)

        if n == 0:  norm_u_accum[n] = norm_u
        else:   norm_u_accum[n] = norm_u_accum[n - 1] + norm_u

        return eu_array, eudiv_array, var_strain_e_u, var_div_e_u, norm_u, norm_u_accum

    # Console output of error and majorant values
    def output_optimization_results(self, iter, maj, m_d, m_f, error, beta, norm_v):
        i_eff = sqrt(maj / error)

        #print('maj opt. # %d, param beta = %.4f: e^2 = %8.2e, maj^2 = %8.2e, m^2_d = %8.2e, m^2_f = %8.2e, i_eff = %.4f' \
        #    % (iter, beta, error, maj, m_d, m_f, i_eff))
        print('maj opt. # %d, param beta = %.4f: e^2 = %8.2e, maj^2 = %8.2e, m^2_d = %8.2e, m^2_f = %8.2e, i_eff = %.4f' \
              % (iter, beta, error / norm_v, maj / norm_v, m_d / norm_v, m_f / norm_v, i_eff))
    # Console output of majorant's values
    def output_maj_optimization_results(self, iter, maj, m_d, m_f, beta):
        print('maj opt. # %d, param beta = %.4f: maj^2 = %8.2e, m^2_d = %8.2e, m^2_f = %8.2e' \
            % (iter, beta, maj, m_d, m_f))
        # print('maj opt. # %d: e^2 = %8.2e, maj^2 = %8.2e, m^2_d = %8.2e, m^2_f = %8.2e, i_eff = %.4f\n' \
        #      % (iter, error, maj, m_d, m_f, i_eff))

    def calculate_majorant_mu_opt(self, u, y, beta, C_FD, f, A, min_eig_A, react):

        # Define optimal parameters
        mu_opt = C_FD ** 2 * (1.0 + beta) * react / (beta * min_eig_A + C_FD ** 2 * (1.0 + beta) * react)

        # Define residuals
        r_d = y - A * grad(u)
        r_f = div(y) + f - react * u

        # Define variational forms
        var_m_d = inner(inv(A) * r_d, r_d)
        var_m_f = inner(r_f, r_f)

        # Define majorant components
        m_d = assemble(var_m_d * dx)
        m_f = assemble(var_m_f * dx)
        m_f_w_opt = assemble(mu_opt / react * var_m_f * dx)  # for calculating majorant
        m_f_one_minus_mu_opt = assemble((1 - mu_opt) ** 2 * var_m_f * dx)  # fo calculating beta_opt

        # Calculate the optimal value for beta parameter
        # beta = C_FD * sqrt(m_f_one_minus_mu_opt / m_d / min_eig_A)
        beta = C_FD * sqrt(m_f / m_d / min_eig_A)

        # Calculate majorant based on the parameter value
        if m_f_w_opt <= DOLFIN_EPS:
            maj = m_d
        else:
            if m_d <= DOLFIN_EPS:
                maj = m_f_w_opt
            else:
                # maj = (1.0 + beta) * m_d + m_f_w_opt
                maj = (1.0 + beta) * m_d + (1.0 + 1 / beta) * C_FD ** 2 * m_f

        return maj, m_d, m_f_one_minus_mu_opt, beta, var_m_d, var_m_f

    def calculate_majorant_bar_mu_opt(self, u, y, beta, C_FD, f_bar, A, min_eig_A, react):

        # Define optimal parameters
        mu_opt = C_FD ** 2 * (1.0 + beta) * react / (beta * min_eig_A + C_FD ** 2 * (1.0 + beta) * react)

        # Define residuals
        r_d = y - A * grad(u)
        r_f = div(y) + f_bar

        # Define variational forms
        #print("inv(A)" , inv(A))

        #var_m_d = abs(inner(inv(A) * r_d, r_d))
        var_m_d = (inner(inv(A) * r_d, r_d))
        var_m_f = inner(r_f, r_f)

        # Define majorant components
        m_d = assemble(var_m_d * dx)
        m_f = assemble(var_m_f * dx)
        #m_f_w_opt_ = assemble(mu_opt**2 / react * var_m_f * dx)  # for calculating majorant
        m_f_w_opt = assemble(mu_opt / react * var_m_f * dx)  # for calculating majorant
        m_f_one_minus_mu_opt = assemble((1 - mu_opt) ** 2 * var_m_f * dx)  # fo calculating beta_opt

        #print("m_f_one_minus_mu_opt = ", m_f_one_minus_mu_opt)
        #print("m_d = ", m_d)
        #print("min_eig_A = ", min_eig_A)
        #input("Press Enter")

        # Calculate majorant based on the parameter value
        if m_f_w_opt <= DOLFIN_EPS:
            maj = m_d
        else:
            if m_d <= DOLFIN_EPS:
                maj = m_f_w_opt
            else:
                maj = (1.0 + beta) * m_d + m_f_w_opt
                #maj_ = (1.0 + beta) * m_d + (1.0 + 1 / beta) * C_FD ** 2 * m_f

        beta = C_FD * sqrt(m_f_one_minus_mu_opt / m_d / min_eig_A)
        # Calculate the optimal value for beta parameter
        #beta_ = C_FD * sqrt(m_f / m_d / min_eig_A)

        return maj, m_d, m_f_w_opt, beta, var_m_d, var_m_f

    # Calculate optimal (wrt to big changes in the reaction fucntion) majorant
    def calculate_majorant(self, u, y, f, mu, lmbda, C_FD, beta):
        """
        :param u: approximate solution
        :param y: flux
        :param beta: parameter minimizing majorant
        :param f_bar: right-hand side function (modified RHS)
        :param A: diffusion matrix
        :param min_eig_A: minimal eigenvalue of diffusion matrix
        :param mesh: mesh
        :param dim: geometrical problem dimension
        :param C_FD: Friedrichs constant
        :return maj: majorant value
                m_d, m_f_one_minus_mu_opt: majorant components
                beta: optimal parameter
                var_m_d, var_m_f_w_opt: variational form of the majorant (to use them later in construction of error estimators)
        """
        # Define residuals
        epsilon = lambda v: 0.5 * (nabla_grad(v) + nabla_grad(v).T)
        sigma = lambda v: 2 * mu * epsilon(v) + lmbda * div(v) * Identity(self.gdim)
        inv_sigma = lambda y: 1 / (2 * mu) * (y - lmbda / (3 * lmbda + 2 * mu) * tr(y) * Identity(self.gdim))

        r_dL = y - sigma(u)
        r_dinvL = inv_sigma(y) - epsilon(u)
        r_f = div(y) + f
        #r_f = nabla_div(y) + f

        # Define variational forms
        # var_m_d = tr(inner(r_dinvL, r_dL))
        var_m_d = tr(inner(inv_sigma(r_dL), r_dL))
        # var_m_d = tr(inner(sigma(u), epsilon(u)))
        #          + 1 / (2 * mu) * (tr(inner(y, y)) - lmbda / (3 * lmbda + 2 * mu) * inner(tr(y), tr(y))) \ # + tr(inner(sigma_inv(y), y))
        #          - 2*tr(inner(epsilon(u), y)))
        var_m_f = inner(r_f, r_f)

        # Define majorant components
        m_d = (assemble(var_m_d * dx))
        m_f = assemble(var_m_f * dx)  # for calculating majorant

        beta = C_FD * sqrt(m_f / m_d)

        # Calculate majorant based on the parameter value
        if m_f <= DOLFIN_EPS:
            maj = m_d
        else:
            if m_d <= DOLFIN_EPS:
                maj = m_f
            else:
                maj = (1.0 + beta) * m_d + (1.0 + 1.0 / beta) * (C_FD ** 2) * m_f
        # Calculate the optimal value for beta parameter
        # beta = C_FD * sqrt(m_f / m_d)

        return maj, m_d, m_f, beta, var_m_d, var_m_f

    def get_matrices_of_optimization_problem(self, H_Div, v, f, mu, lmbda):
        """
        :param H_div: funtional space for the flux function
        :param v: approximate solution
        :param f_bar: right-hand side function (modified RHS)
        :param A: diffusion matrix
        :param invA: diffusion matrix
        :return divdiv, PhiPhi, RhsdivPhi, RhsPhi: matrixes for the optimal reconstruction of the flux
        """
        # Define variational problem
        y = TrialFunction(H_Div)  # tensor-field
        z = TestFunction(H_Div)  # tensor-field

        epsilon = lambda v: 0.5 * (nabla_grad(v) + nabla_grad(v).T)
        inv_sigma = lambda y: 1 / (2 * mu) * (y - lmbda / (3 * lmbda + 2 * mu) * tr(y) * Identity(self.gdim))

        # Define system of linear equation to find the majorant
        # Phi_i, i = 1, ..., d are vector-valued basis of H_div
        # \int_\Omega div(phi_i) div(phi_j) \dx
        divPhidivPhi = assemble(inner(div(y), div(z)) * dx)
        #divPhidivPhi = assemble(inner(nabla_div(y), nabla_div(z)) * dx(domain=self.mesh))
        # \int_\Omega L^{-1} phi_i \cdot phi_j \dx
        PhiPhi = assemble(tr(inner(inv_sigma(y), z)) * dx)
        # Define vectors included into the RHS
        RhsdivPhi = assemble(inner(-f, div(z)) * dx)
        #RhsdivPhi = assemble(inner(-f, nabla_div(z)) * dx(domain=self.mesh))
        RhsPhi = assemble(tr(inner(epsilon(v), z)) * dx)

        return divPhidivPhi, PhiPhi, RhsdivPhi, RhsPhi

    def get_matrices_of_optimization_problem_nobar(self, H_div, v, f, A, react):
        """
        :param H_div: funtional space for the flux function
        :param v: approximate solution
        :param f_bar: right-hand side function (modified RHS)
        :param A: diffusion matrix
        :param invA: diffusion matrix
        :return divdiv, PhiPhi, RhsdivPhi, RhsPhi: matrixes for the optimal reconstruction of the flux
        """
        # Define variational problem
        y = TrialFunction(H_div)
        q = TestFunction(H_div)

        # Define system of linear equation to find the majorant
        # Phi_i, i = 1, ..., d are vector-valued basis of H_div
        # \int_\Omega div(phi_i) div(phi_j) \dx
        divPhidivPhi = assemble(inner(nabla_div(y), nabla_div(q)) * dx(domain=self.mesh))
        # \int_\Omega A^{-1} phi_i \cdot phi_j \dx
        PhiPhi = assemble(inner(inv(A) * y, q) * dx(domain=self.mesh))
        # Define vectors included into the RHS
        RhsdivPhi = assemble(inner(-(f - react*v), nabla_div(q)) * dx(domain=self.mesh))
        RhsPhi = assemble(inner(nabla_grad(v), q) * dx(domain=self.mesh))

        return divPhidivPhi, PhiPhi, RhsdivPhi, RhsPhi

    def get_matrices_of_optimization_problem_bar(self, H_div, v, f_bar, A):
        """
        :param H_div: funtional space for the flux function
        :param v: approximate solution
        :param f_bar: right-hand side function (modified RHS)
        :param A: diffusion matrix
        :param invA: diffusion matrix
        :return divdiv, PhiPhi, RhsdivPhi, RhsPhi: matrixes for the optimal reconstruction of the flux
        """
        # Define variational problem
        y = TrialFunction(H_div)
        q = TestFunction(H_div)

        # Define system of linear equation to find the majorant
        # Phi_i, i = 1, ..., d are vector-valued basis of H_div
        # \int_\Omega div(phi_i) div(phi_j) \dx
        divPhidivPhi = assemble(inner(div(y), div(q)) * dx) #divPhidivPhi = assemble(inner(nabla_div(y), nabla_div(q)) * dx)
        # \int_\Omega A^{-1} phi_i \cdot phi_j \dx
        PhiPhi = assemble(inner(inv(A) * y, q) * dx)
        # Define vectors included into the RHS
        RhsdivPhi = assemble(inner(-f_bar, div(q)) * dx)    #RhsdivPhi = assemble(inner(-f_bar, nabla_div(q)) * dx)
        RhsPhi = assemble(inner(grad(v), q) * dx)           #RhsPhi = assemble(inner(nabla_grad(v), q) * dx)

        return divPhidivPhi, PhiPhi, RhsdivPhi, RhsPhi

    def majorant_pressure(self, maj_ep_array, maj_epl2_array, ep_array,
                          k, p, f, react,
                          domain_params, test_params,
                          func_spaces, funcs, norm_p):

        C_FD = domain_params["C_FD"]

        # Initial flux guess
        #y = project(funcs["kappa"] * nabla_grad(p), func_spaces["H_div"])
        #y = project(funcs["kappa"] * grad(p), func_spaces["H_div"])
        #print("kappa = ", funcs["kappa"])
        y = project(funcs["kappa"] * grad(interpolate(funcs["p_e"], self.Qe)), func_spaces["H_div"]) # exact flux for debugging
        # Auxiliary Young's parameter
        beta = 1.0

        # f_bar is new RHS
        f_bar = f - react * p
        '''
        maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = \
            self.calculate_majorant_mu_opt(p, y, beta, C_FD,
                                           f, funcs["kappa"], funcs["min_eig_kappa"],
                                           react)
        '''
        maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = \
            self.calculate_majorant_bar_mu_opt(p, y, beta, C_FD,
                                               f_bar, funcs["kappa"], funcs["min_eig_kappa"],
                                               react)

        print('%------------------------------------------------------------------------------------%')
        print("% Majorant for e_p")
        print('%------------------------------------------------------------------------------------%')

        self.output_optimization_results(0, maj, m_d, m_f, ep_array[k - 1], beta, norm_p)

        if test_params["majorant_optimisation"]:
            # ----------------------------------------------------------------------------#
            # Optimization algorithm
            # ----------------------------------------------------------------------------#
            '''
            S, K, z, g = self.get_matrices_of_optimization_problem_nobar(func_spaces["H_div"], p,
                                                                       f_bar, funcs["kappa"], react)
            '''
            S, K, z, g = self.get_matrices_of_optimization_problem_bar(func_spaces["H_div"], p,
                                                                       f_bar, funcs["kappa"])
            y = Function(func_spaces["H_div"])
            Y = y.vector()

            # Execute iterative process to optimize majorant with respect to beta and flux
            for i in range(0, test_params["majorant_optimization_iterations_number"]):
                # Define the optimal system for the flux reconstruction
                yMatrix = C_FD ** 2 / funcs["min_eig_kappa"] * S + beta * K
                yRhs = C_FD ** 2 / funcs["min_eig_kappa"] * z + beta * g
                solve(yMatrix, Y, yRhs)
                y.vector = Y
                maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = \
                    self.calculate_majorant_bar_mu_opt(p, y, beta, C_FD,
                                                       f_bar, funcs["kappa"], funcs["min_eig_kappa"],
                                                       react)

                self.output_optimization_results(i + 1, maj, m_d, m_f, ep_array[k - 1], beta, norm_p)

        maj_ep_array[k - 1] = maj
        maj_epl2_array[k - 1] = (C_FD ** (-2) * funcs["min_eig_kappa"] + react) ** (-1) * maj

        return maj_ep_array, maj_epl2_array, var_m_d

    def majorant_pressure_part(self,
                               n, maj_d, maj_eq,
                               p, f,
                               domain_params, test_params,
                               func_spaces, funcs):

        C_FD = domain_params["C_FD"]

        # Initial flux guess
        #y = project(funcs["kappa"] * nabla_grad(p), func_spaces["H_div"])
        y = project(funcs["kappa"] * grad(p), func_spaces["H_div"])
        #y = project(funcs["kappa"] * grad(interpolate(funcs["p_e"], self.Qe)), func_spaces["H_div"]) # exact flux for debugging
        # Auxiliary Young's parameter
        beta = 1.0

        # f_bar is new RHS
        maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = \
            self.calculate_majorant_bar_mu_opt(p, y, beta, C_FD, f, funcs["kappa"], funcs["min_eig_kappa"], funcs["react"])
        #'''
        print('%------------------------------------------------------------------------------------%')
        print("% Minimising M^2_h(p)")

        self.output_maj_optimization_results(0, maj, m_d, m_f, beta)

        if test_params["majorant_optimisation"]:
            # ----------------------------------------------------------------------------#
            # Optimization algorithm
            # ----------------------------------------------------------------------------#
            S, K, z, g = self.get_matrices_of_optimization_problem_bar(func_spaces["H_div"], p,
                                                                       f, funcs["kappa"])
            y = Function(func_spaces["H_div"])
            Y = y.vector()

            # Execute iterative process to optimize majorant with respect to beta and flux
            for i in range(0, test_params["majorant_optimization_iterations_number"]):
                # Define the optimal system for the flux reconstruction
                yMatrix = C_FD ** 2 / funcs["min_eig_kappa"] * S + beta * K
                yRhs = C_FD ** 2 / funcs["min_eig_kappa"] * z + beta * g
                solve(yMatrix, Y, yRhs)
                y.vector = Y
                maj, m_d, m_f, beta, var_m_d, var_m_f_w_opt = \
                    self.calculate_majorant_bar_mu_opt(p, y, beta, C_FD, f, funcs["kappa"], funcs["min_eig_kappa"], funcs["react"])

                self.output_maj_optimization_results(i + 1, maj, m_d, m_f, beta)

        maj_d[n - 1] = m_d
        maj_eq[n - 1] = m_f

        return maj_d, maj_eq, var_m_d

    def majorant_displacement(self, maj_eu_array, eu_array,
                              k, u, F, mu, lmbda,
                              domain_params, test_params,
                              func_spaces,
                              u_e, norm_u):

        C_K = math.sqrt(2) * domain_params["C_FD"] / sqrt(2 * mu)
        beta = 1.0

        # Initial flux guess
        # y = project(sigma(u, mu, lmbda, self.gdim), func_spaces["H_Div"])
        epsilon = lambda v: 0.5 * (nabla_grad(v) + nabla_grad(v).T)
        sigma = lambda v: 2 * mu * epsilon(v) + lmbda * div(v) * Identity(self.gdim)
        y = project(sigma(interpolate(u_e, self.VVe)), func_spaces["H_Div"])

        # beta is the reaction fucntion
        maj, m_d, m_f, beta, var_m_d, var_m_f = \
            self.calculate_majorant(u, y, F, mu, lmbda, C_K, beta)

        print('%------------------------------------------------------------------------------------%')
        print("% Majorant for e_u")
        print('%------------------------------------------------------------------------------------%')

        self.output_optimization_results(0, maj, m_d, m_f, eu_array[k - 1], beta, norm_u)

        if test_params["majorant_optimisation"]:
            # ----------------------------------------------------------------------------#
            # Optimization algorithm
            # ----------------------------------------------------------------------------#
            S, K, z, g = self.get_matrices_of_optimization_problem(func_spaces["H_Div"], u, F, mu, lmbda)

            y = Function(func_spaces["H_Div"])
            Y = y.vector()

            old_maj = DOLFIN_EPS

            # Execute iterative process to optimize majorant with respect to beta and flux
            for i in range(0, 0):
                yMatrix = C_K ** 2 * S + beta * K
                yRhs = C_K ** 2 * z + beta * g

                solve(yMatrix, Y, yRhs)
                y.vector = Y

                # Calculate majorant
                maj, m_d, m_f, beta, var_m_d, var_m_f = \
                    self.calculate_majorant(u, y, F, mu, lmbda, C_K, beta)

                self.output_optimization_results(i + 1, maj, m_d, m_f, eu_array[k - 1], beta, norm_u)
                # print("majorant improvement: %10.4e\n" % (abs(old_maj - maj) / old_maj))
                old_maj = maj

        maj_eu_array[k - 1] = maj

        return maj_eu_array

    def majorant_displacement_part(self, n, maj_d, maj_eq,
                                  u, F, mu, lmbda,
                                  domain_params, test_params,
                                  func_spaces, funcs):

        C_K = math.sqrt(2) * domain_params["C_FD"] / sqrt(2 * mu)
        beta = 1.0

        # Initial flux guess
        # y = project(sigma(u, mu, lmbda, self.gdim), func_spaces["H_Div"])
        epsilon = lambda v: 0.5 * (nabla_grad(v) + nabla_grad(v).T)
        sigma = lambda v: 2 * mu * epsilon(v) + lmbda * div(v) * Identity(self.gdim)
        #y = project(sigma(interpolate(u, self.VVe)), func_spaces["H_Div"])
        y = project(sigma(interpolate(funcs["u_e"], self.VVe)), func_spaces["H_Div"])

        # beta is the reaction fucntion
        maj, m_d, m_f, beta, var_m_d, var_m_f = \
            self.calculate_majorant(u, y, F, mu, lmbda, C_K, beta)

        print('%------------------------------------------------------------------------------------%')
        print("% Minimisation of M^2(u)")
        print('%------------------------------------------------------------------------------------%')

        self.output_maj_optimization_results(0, maj, m_d, m_f, beta)

        if test_params["majorant_optimisation"]:
            # ----------------------------------------------------------------------------#
            # Optimization algorithm
            # ----------------------------------------------------------------------------#
            S, K, z, g = self.get_matrices_of_optimization_problem(func_spaces["H_Div"], u, F, mu, lmbda)

            y = Function(func_spaces["H_Div"])
            Y = y.vector()

            # Execute iterative process to optimize majorant with respect to beta and flux
            for i in range(0, 0):
                yMatrix = C_K ** 2 * S + beta * K
                yRhs = C_K ** 2 * z + beta * g

                solve(yMatrix, Y, yRhs)
                y.vector = Y

                # Calculate majorant
                maj, m_d, m_f, beta, var_m_d, var_m_f = \
                    self.calculate_majorant(u, y, F, mu, lmbda, C_K, beta)

                self.output_optimization_results(i + 1, maj, m_d, m_f, eu_array[k - 1], beta)
                # print("majorant improvement: %10.4e\n" % (abs(old_maj - maj) / old_maj))
                old_maj = maj

        maj_d[k - 1] = m_d
        maj_eq[k - 1] = m_f
        return maj_d, maj_eq, var_m_d


# Class to solve the Biot problem
class BiotSolvers():

    # Init the BiotSolver with the mesh
    def __init__(self, mesh, material_params, domain_params, test_params):
        self.mesh = mesh
        self.material_params = material_params
        self.domain_params = domain_params
        self.test_params = test_params

    '''
    def get_last_iter(self, n, i,
                      ep_enrg, ep_l2, eu_enrg, eu_div, maj_ph, maj_phl2, maj_uh, maj_uhdiv, maj_pi, maj_ui,
                      ep_enrg_it, ep_l2_it, eu_enrg_it, eu_div_it, maj_ph_it, maj_phl2_it, maj_uh_it, maj_uhdiv_it,
                      maj_pi_it, maj_ui_it, norm_p, norm_u, norm_p_accum, norm_u_accum):
        if n == 0 :
            norm_p_accum[n] = norm_p
            norm_u_accum[n] = norm_u

            ep_enrg[n] = ep_enrg_it[i] / norm_p
            ep_l2[n]   = ep_l2_it[i] / norm_p
            eu_enrg[n] = eu_enrg_it[i] / norm_u
            eu_div[n]  = eu_div_it[i] / norm_u

            maj_ph[n] = maj_ph_it[i] / norm_p
            maj_phl2[n] = maj_phl2_it[i] / norm_p
            maj_uh[n] = maj_uh_it[i] / norm_u
            maj_uhdiv[n] = maj_uhdiv_it[i] / norm_u
            maj_pi[n] = maj_pi_it[i] / norm_p
            maj_ui[n] = maj_ui_it[i] / norm_u
        else:
            norm_p_accum[n] = norm_p_accum[n - 1] + norm_p
            norm_u_accum[n] = norm_u_accum[n - 1] + norm_u

            ep_enrg[n] = (ep_enrg[n-1] + ep_enrg_it[i]) / norm_p_accum[n]
            ep_l2[n] = (ep_l2[n-1] + ep_l2_it[i]) / norm_p_accum[n]
            eu_enrg[n] = (eu_enrg[n-1] + eu_enrg_it[i]) / norm_u_accum[n]
            eu_div[n] = (eu_div[n-1] + eu_div_it[i]) / norm_u_accum[n]

            maj_ph[n] = (maj_ph[n-1] + maj_ph_it[i]) / norm_p_accum[n]
            maj_phl2[n] = (maj_phl2[n-1] + maj_phl2_it[i]) / norm_p_accum[n]
            maj_uh[n] = (maj_uh[n-1] + maj_uh_it[i]) / norm_u_accum[n]
            maj_uhdiv[n] = (maj_uhdiv[n-1] + maj_uhdiv_it[i]) / norm_u_accum[n]

            maj_pi[n] = (maj_pi[n-1] + maj_pi_it[i]) / norm_p_accum[n]
            maj_ui[n] = (maj_ui[n-1] + maj_ui_it[i]) / norm_u_accum[n]

        return ep_enrg, ep_l2, eu_enrg, eu_div, maj_ph, maj_phl2, maj_uh, maj_uhdiv, maj_pi, maj_ui, norm_p_accum, norm_u_accum
    '''

    def get_last_iter(self, n, i,
                      ep_enrg, ep_l2, eu_enrg, eu_div, maj_ph, maj_phl2, maj_uh, maj_uhdiv, maj_pi, maj_ui,
                      ep_enrg_it, ep_l2_it, eu_enrg_it, eu_div_it,
                      maj_ph_it, maj_phl2_it, maj_uh_it, maj_uhdiv_it,
                      maj_pi_it, maj_ui_it,
                      maj_tilde_pi_it, maj_tilde_ui_it,
                      maj_pi_it_aux, maj_ui_it_aux,
                      e, e_p, e_u, maj_p, maj_u, maj,
                      norm_p, norm_u,
                      norm_p_accum, norm_u_accum):
        if n == 0:
            norm_p_accum[n] = norm_p
            norm_u_accum[n] = norm_u

            ep_enrg[n] = ep_enrg_it[i]
            ep_l2[n] = ep_l2_it[i]
            eu_enrg[n] = eu_enrg_it[i]
            eu_div[n] = eu_div_it[i]

            maj_ph[n] = maj_ph_it[i]
            maj_phl2[n] = maj_phl2_it[i]
            maj_uh[n] = maj_uh_it[i]
            maj_uhdiv[n] = maj_uhdiv_it[i]

            maj_pi[n] = maj_pi_it[i]
            maj_ui[n] = maj_ui_it[i]

            # Correction of the maj_pi and maj_ui
            maj_pi[n] = maj_pi_it_aux
            maj_ui[n] = maj_ui_it_aux

            # Correction of the maj_pi and maj_ui by maj_tilde_pi and maj_tilde_ui
            # maj_pi[n] = maj_tilde_pi_it[i]
            # maj_ui[n] = maj_tilde_ui_it[i]

            # Using auxiliary majorant
            e_p[n] = ep_enrg_it[i]
            e_u[n] = eu_enrg_it[i]
            e[n] = ep_enrg_it[i] + eu_enrg_it[i]

            maj_p[n] = 2 * (maj_pi_it_aux + maj_ph_it[i])
            maj_u[n] = 2 * (maj_ui_it_aux + maj_uh_it[i])
            maj[n] = 2 * (maj_pi_it_aux + maj_ph_it[i] + maj_ui_it_aux + maj_uh_it[i])

        else:

            #print("norm_p_accum[n - 1]", norm_p_accum[n - 1])
            #print("norm_p_accum[n - 1]", norm_p_accum[n - 1])
            #print("norm_p_accum[n]", norm_p_accum[n])
            #print("norm_u_accum[n]", norm_u_accum[n])

            norm_p_accum[n] = norm_p_accum[n - 1] + norm_p
            norm_u_accum[n] = norm_u_accum[n - 1] + norm_u

            ep_enrg[n] = ep_enrg[n - 1] + ep_enrg_it[i]
            ep_l2[n] = ep_l2[n - 1] + ep_l2_it[i]
            eu_enrg[n] = eu_enrg[n - 1] + eu_enrg_it[i]
            eu_div[n] = eu_div[n - 1] + eu_div_it[i]

            maj_ph[n] = maj_ph[n - 1] + maj_ph_it[i]
            maj_phl2[n] = maj_phl2[n - 1] + maj_phl2_it[i]
            maj_uh[n] = maj_uh[n - 1] + maj_uh_it[i]
            maj_uhdiv[n] = maj_uhdiv[n - 1] + maj_uhdiv_it[i]

            maj_pi[n] = maj_pi[n - 1] + maj_pi_it[i]
            maj_ui[n] = maj_ui[n - 1] + maj_ui_it[i]

            # Using auxiliary majorant
            e_p[n] = e_p[n - 1] + ep_enrg_it[i]
            e_u[n] = e_u[n - 1] + eu_enrg_it[i]
            e[n] = e[n - 1] + ep_enrg_it[i] + eu_enrg_it[i]

            maj_p[n] = maj_p[n - 1] + 2 * (maj_pi_it_aux + maj_ph_it[i])
            maj_u[n] = maj_u[n - 1] + 2 * (maj_ui_it_aux + maj_uh_it[i])
            maj[n] = maj[n - 1] + 2 * (maj_pi_it_aux + maj_ph_it[i] + maj_ui_it_aux + maj_uh_it[i])

        '''
        if n == 0 :
            norm_p_accum[n] = norm_p
            norm_u_accum[n] = norm_u

            ep_enrg[n] = ep_enrg_it[i] / norm_p
            ep_l2[n] = ep_l2_it[i] / norm_p
            eu_enrg[n] = eu_enrg_it[i] / norm_u
            eu_div[n] = eu_div_it[i] / norm_u

            maj_ph[n] = maj_ph_it[i] / norm_p
            maj_phl2[n] = maj_phl2_it[i] / norm_p
            maj_uh[n] = maj_uh_it[i] / norm_u
            maj_uhdiv[n] = maj_uhdiv_it[i] / norm_u


            maj_pi[n] = maj_pi_it[i] / norm_p
            maj_ui[n] = maj_ui_it[i] / norm_u

            # Correction of the maj_pi and maj_ui
            maj_pi[n] = maj_pi_it_aux / norm_p_accum[n]
            maj_ui[n] = maj_ui_it_aux / norm_u_accum[n]

            # Using auxiliary majorant
            e_p[n] = ep_enrg_it[i] / norm_p_accum[n]
            e_u[n] = eu_enrg_it[i] / norm_u_accum[n]
            e[n] = (ep_enrg_it[i] + eu_enrg_it[i]) / max([norm_p_accum[n], norm_u_accum[n]])

            maj_p[n] = 2 * (maj_pi_it_aux + maj_ph_it[i]) / norm_p_accum[n]
            maj_u[n] = 2 * (maj_ui_it_aux + maj_uh_it[i]) / norm_u_accum[n]
            maj[n] = 2 * (maj_pi_it_aux + maj_ph_it[i] + maj_ui_it_aux + maj_uh_it[i]) / max([norm_p_accum[n], norm_u_accum[n]])


        else:

            print("norm_p_accum[n - 1]", norm_p_accum[n - 1])
            print("norm_p_accum[n - 1]", norm_p_accum[n - 1])
            print("norm_p_accum[n]", norm_p_accum[n])
            print("norm_u_accum[n]", norm_u_accum[n])

            norm_p_accum[n] = norm_p_accum[n - 1] + norm_p
            norm_u_accum[n] = norm_u_accum[n - 1] + norm_u

            ep_enrg[n] = (ep_enrg[n-1] + ep_enrg_it[i]) / norm_p_accum[n]
            ep_l2[n] = (ep_l2[n-1] + ep_l2_it[i]) / norm_p_accum[n]
            eu_enrg[n] = (eu_enrg[n-1] + eu_enrg_it[i]) / norm_u_accum[n]
            eu_div[n] = (eu_div[n-1] + eu_div_it[i]) / norm_u_accum[n]

            maj_ph[n] = (maj_ph[n-1] + maj_ph_it[i]) / norm_p_accum[n]
            maj_phl2[n] = (maj_phl2[n-1] + maj_phl2_it[i]) / norm_p_accum[n]
            maj_uh[n] = (maj_uh[n-1] + maj_uh_it[i]) / norm_u_accum[n]
            maj_uhdiv[n] = (maj_uhdiv[n-1] + maj_uhdiv_it[i]) / norm_u_accum[n]

            maj_pi[n] = (maj_pi[n-1] + maj_pi_it[i]) / norm_p_accum[n]
            maj_ui[n] = (maj_ui[n-1] + maj_ui_it[i]) / norm_u_accum[n]

            # Using auxiliary majorant
            e_p[n] = (e_p[n-1] + ep_enrg_it[i]) / norm_p_accum[n]
            e_u[n] = (e_u[n-1] + eu_enrg_it[i]) / norm_u_accum[n]
            e[n] = (e[n-1] + ep_enrg_it[i] + eu_enrg_it[i]) / max([norm_p_accum[n], norm_u_accum[n]])

            maj_p[n] = (maj_p[n-1] + 2 * (maj_pi_it_aux + maj_ph_it[i])) / norm_p_accum[n]
            maj_u[n] = (maj_u[n-1] + 2 * (maj_ui_it_aux + maj_uh_it[i])) / norm_u_accum[n]
            maj[n] = (maj[n-1] + 2 * (maj_pi_it_aux + maj_ph_it[i] + maj_ui_it_aux + maj_uh_it[i])) / max([norm_p_accum[n], norm_u_accum[n]])

        print("ep_enrg_it[i] + eu_enrg_it[i] = ", ep_enrg_it[i] + eu_enrg_it[i])
        print("2 * (maj_pi_it_aux + maj_ph_it[i] + maj_ui_it_aux + maj_uh_it[i]) = ", 2 * (maj_pi_it_aux + maj_ph_it[i] + maj_ui_it_aux + maj_uh_it[i]))

        input("Press")
        '''
        return ep_enrg, ep_l2, eu_enrg, eu_div, \
               maj_ph, maj_phl2, maj_uh, maj_uhdiv, maj_pi, maj_ui, \
               e, e_p, e_u, \
               maj_p, maj_u, maj, \
               norm_p_accum, norm_u_accum

    # Function with iterative coupling procedure
    def iterative_coupling(self, n, u, p, v, theta, bc_p, bc_u,
                           error_control,
                           ep_enrg, ep_l2, eu_enrg, eu_div, e_p, e_u, e,
                           maj_ph, maj_phl2, maj_uh, maj_uhdiv,
                           maj_pi, maj_ui,
                           maj_p, maj_u, maj,
                           alpha, mu, lmbda, L, q, ds,
                           u_i, u_i1, p_i, p_i1,
                           g_tilde, funcs, func_spaces,
                           test_params,
                           norm_p_accum, norm_u_accum):

        gdim = domain_params['gdim']
        C_FD = domain_params["C_FD"]

        ep_enrg_it = error_control.get_error_array(test_params["iter_num"])
        eu_enrg_it = error_control.get_error_array(test_params["iter_num"])
        ep_l2_it = error_control.get_error_array(test_params["iter_num"])
        eu_div_it = error_control.get_error_array(test_params["iter_num"])

        maj_ph_it = error_control.get_error_array(test_params["iter_num"])
        maj_phl2_it = error_control.get_error_array(test_params["iter_num"])
        maj_uh_it = error_control.get_error_array(test_params["iter_num"])
        maj_uhdiv_it = error_control.get_error_array(test_params["iter_num"])
        maj_pi_it = error_control.get_error_array(test_params["iter_num"])
        maj_ui_it = error_control.get_error_array(test_params["iter_num"])

        maj_tilde_pi_it = error_control.get_error_array(test_params["iter_num"])
        maj_tilde_ui_it = error_control.get_error_array(test_params["iter_num"])

        maj_pi_it_mult = error_control.get_error_array(test_params["iter_num"])
        maj_ui_it_mult = error_control.get_error_array(test_params["iter_num"])
        q_it = error_control.get_error_array(test_params["iter_num"])
        maj_base_it = error_control.get_error_array(test_params["iter_num"])

        incr_sigma_norm_it = error_control.get_error_array(test_params["iter_num"])
        incr_p_norm_it = error_control.get_error_array(test_params["iter_num"])
        incr_u_norm_it = error_control.get_error_array(test_params["iter_num"])

        gamma = math.sqrt(2 * L)

        u_0h = u_i
        p_0h = p_i
        
        u_0 = interpolate(funcs['u_0'], func_spaces['V_exact'])
        p_0 = interpolate(funcs['p_0'], func_spaces['Q_exact'])

        eta_0h = alpha / gamma * div(u_0h) - L / gamma * p_0h
        eta_0 = alpha / gamma * div(u_0) - L / gamma * p_0
        diff_eta = eta_0 - eta_0
        norm_diff_eta_0 = assemble(inner(diff_eta, diff_eta) * dx)

        print("norm_diff_eta_0 = ", norm_diff_eta_0)


        for i in range(1, test_params["iter_num"] + 1):

            if test_params["full_documentation"]:
                print('%---------------------------------------------------------------------------------%')
                print('% Iteration # ', i)

            epsilon = lambda v: 0.5 * (nabla_grad(v) + nabla_grad(v).T)
            sigma = lambda v: 2 * mu * epsilon(v) + lmbda * div(v) * Identity(gdim)

            # ----------------------------------------------------------------------------------------------#
            # step 1: Constructing p (p_{i+1}) using u_i
            # ----------------------------------------------------------------------------------------------#
            a_p_var = (funcs["react"] * inner(p, theta) + inner(funcs["kappa"] * grad(p), grad(theta))) * dx
            G_i = L * p_i - alpha * div(u_i) + g_tilde
            #G_i = L * p_i - alpha * tr(epsilon(u_i)) + g_tilde
            l_p_var = inner(G_i, theta) * dx

            # Assemble the system to solve flow equation for pressure
            A_p, b_p = assemble_system(a_p_var, l_p_var, bc_p)
            solve(A_p, p_i1.vector(), b_p)
            # solve(A_p, p_k1.vector(), b_p, "gmres", "hypre_amg")

            # Calculate the error in the pressure term
            ep_enrg_it, ep_l2_it, var_grad_ep, var_ep, norm_p, norm_p_accum \
                = error_control.error_pressure(ep_enrg_it, ep_l2_it, i, p_i1, norm_p_accum, n, funcs)

            # ----------------------------------------------------------------------------------------------#
            # Constructin u (u_{i+1}) using p_{i+1}
            # ----------------------------------------------------------------------------------------------#
            F_i = funcs["f"] - alpha * grad(p_i1)
            # tr here to map ComponentTensor(2, 2) to Trace that can be integrated
            a_u_var = tr(inner(sigma(u), epsilon(v))) * dx
            l_u_var = inner(F_i, v) * dx + inner(funcs["t_N"], v) * ds(4)

            # Assemble the system to solve mechanics equation for displacement
            A_u, b_u = assemble_system(a_u_var, l_u_var, bc_u)
            solve(A_u, u_i1.vector(), b_u)

            # TODO: add here print out of inc_p = || p^{k+1} - p^k ||_l2 and inc_u = || u^{k+1} - u^k ||_l2
            # TODO: add the exiting criterion based on the inc_p, inc_u <= 1e-8

            # Calculate the error in the displacement term
            eu_enrg_it, eu_div_it, var_strain_eu, var_div_eu, norm_u, norm_u_accum \
                = error_control.error_displacement(eu_enrg_it, eu_div_it, i, u_i1, norm_u_accum, n, funcs, mu, lmbda)


            eta_i1 = project(alpha / gamma * div(u_i1) - L / gamma * p_i1, func_spaces["Q"])
            eta_i = project(alpha / gamma * div(u_i) - L / gamma * p_i, func_spaces["Q"])

            if i == 1:
                incr_sigma_0_norm = assemble((eta_i1 - eta_i) ** 2 * dx)

            pow = test_params["pow"]
            aux_it = test_params["iter_num"] - pow

            if i == aux_it:
                # Make a deep copy of eta_i
                eta_aux = eta_i.copy(deepcopy=True)

            if self.test_params["error_estimates"] == True and ((i >= test_params["iter_num"] - 1 and i <= test_params["iter_num"]) or i == aux_it):

                # Reconstruct pressure majorant
                maj_ph_it, maj_phl2_it, var_md_ep = \
                    error_control.majorant_pressure(maj_ph_it, maj_phl2_it, ep_enrg_it,
                                                    i, p_i1, G_i, funcs["react"],
                                                    self.domain_params, self.test_params,
                                                    func_spaces, funcs, norm_p)
                """
                ep_distr, md_ep_distr = \
                    error_control.error_majorant_distribution_pressure(func_spaces,
                                                                       var_ep, var_md_ep)
                """
                if test_params["full_documentation"]:
                    print("||| e_p |||^2 = %.4e" % (ep_enrg_it[i-1] / norm_p))
                    print("M^{h,2}(p)    = %.4e" % (maj_ph_it[i-1] / norm_p))
                    print("|| e_p ||^2   = %.4e" % (ep_l2_it[i-1] / norm_p))
                    print("M^{h,2}_L2(p) = %.4e\n" % (maj_phl2_it[i-1] / norm_p))

                # Reconstruct pressure majorant
                maj_uh_it = \
                    error_control.majorant_displacement(maj_uh_it, eu_enrg_it,
                                                        i, u_i1, F_i, mu, lmbda,
                                                        self.domain_params, self.test_params,
                                                        func_spaces,
                                                        funcs["u_e"], norm_u)
                nu = 1 / (lmbda)  # Young's parameter, nu >= 1/(2*lmbda), e.g., nu = 2/(3*lmbda), 3/(4*lmbda)
                maj_uh_it[i-1] = 2 * (
                        maj_uh_it[i-1] + lmbda * (nu * alpha) ** 2 / (2 * nu * lmbda - 1) * maj_phl2_it[i-1])
                maj_uhdiv_it[i-1] = maj_uh_it[i-1] / (2 * mu * gdim + lmbda)

                if test_params["full_documentation"]:
                    print("||| e_u |||^2    = %.4e" % (eu_enrg_it[i-1] / norm_u) )
                    print("M^{h,2}(u)       = %.4e" % (maj_uh_it[i-1] / norm_u))
                    print("|| div(e_u) ||^2 = %.4e" % (eu_div_it[i-1] / norm_u))
                    print("M^{h,2}_div(u)   = %.4e\n" % (maj_uhdiv_it[i-1] / norm_u))

                if i == test_params["iter_num"]:

                    incr_sigma_norm = assemble((eta_i1 - eta_i) ** 2 * dx)

                    #incr_p_norm = assemble((p_i1 - p_i) ** 2 * dx)
                    #incr_divu_norm = assemble(div(u_i1 - u_i) ** 2 * dx)

                    #print("\| p_i1 - p_i \|      = %.4e" % (incr_p_norm))
                    #print("\| div(u_i1 - u_i) \| = %.4e\n" % (incr_divu_norm))

                    print("Components of the iterative majorants:")
                    print("--------------------------------------")

                    maj_base = incr_sigma_norm \
                               + lmbda / 2 * (maj_uhdiv_it[i-1] + maj_uhdiv_it[i - 2]) \
                               + L / 4 * (maj_phl2_it[i-1] + maj_phl2_it[i - 2])

                    maj_tilde_base = q**(2 * (test_params["iter_num"] - 1)) * ((q**2 + 1) * incr_sigma_0_norm + norm_diff_eta_0)

                    print("\| eta_i1 - eta_i \|  = %.4e" % (incr_sigma_norm))
                    print("3 * q / (1 - q)**2 * (C_FD ** 2 * beta / min_eig_kappa + 1) = %.4e, " %
                          (3 * q / (1 - q) ** 2 * (C_FD ** 2 * funcs["beta"] / funcs["min_eig_kappa"] + 1)))
                    #print("where")
                    #print("(C_FD ** 2 * beta / min_eig_kappa + 1) = %.4e" %
                    #      ((C_FD ** 2 * funcs["beta"] / funcs["min_eig_kappa"] + 1)))
                    #print("3 * q / (1 - q)**2 = %.4e" %
                    #      (3 * q / (1 - q) ** 2))
                    #print("3 * q**2 / (1 - q)**2 * (1 + gdim * lmbda / (2 * mu)) = %.4e" %
                    #      (3 * q ** 2 / (1 - q) ** 2 * (1 + gdim * lmbda / (2 * mu))))
                    print("maj(i, i-1) = %.4e\n " % maj_base)
                    print("maj_tilde_base = %.4e\n " % maj_tilde_base)

                    maj_pi_it[i-1] = maj_base * 3 * q / (1 - q) ** 2 * (C_FD ** 2 * funcs["beta"] / funcs["min_eig_kappa"] + 1)
                    maj_ui_it[i-1] = maj_base * 3 * q ** 2 / (1 - q) ** 2 * (1 + gdim * lmbda / (2 * mu))

                    maj_tilde_pi_it[i-1] = maj_tilde_base * 3 * q / (1 - q) ** 2 * (
                                C_FD ** 2 * funcs["beta"] / funcs["min_eig_kappa"] + 1)
                    maj_tilde_ui_it[i-1] = maj_tilde_base * 3 * q ** 2 / (1 - q) ** 2 * (1 + gdim * lmbda / (2 * mu))

                    if test_params["full_documentation"]:
                        print("M^{2,i}(p) = %.4e" % (maj_pi_it[i-1] / norm_p))
                        print("M^{2,i}(u) = %.4e" % (maj_ui_it[i-1] / norm_u))

                        print("M_tilde^{2,i}(p) = %.4e" % (maj_tilde_pi_it[i - 1] / norm_p))
                        print("M_tilde^{2,i}(u) = %.4e" % (maj_tilde_ui_it[i - 1] / norm_u))

            # Update the iterations
            u_i.assign(u_i1)
            p_i.assign(p_i1)

        last_it = test_params["iter_num"] - 1

        # Output the convergence of the pressure and displacement
        postprocess.output_errors_and_estimates_wrt_iterations(ep_enrg_it, ep_l2_it, eu_enrg_it, eu_div_it,
                                                               maj_ph_it, maj_phl2_it, maj_uh_it, maj_uhdiv_it,
                                                               maj_pi_it, maj_ui_it,
                                                               test_params["iter_num"], norm_p, norm_u)
        #incr_sigma_norm   = assemble((eta_i1 - eta_i) ** 2 * dx)
        #incr_sigma_norm_0 = assemble((eta_i1 - eta_0) ** 2 * dx)
        #incr_p_norm_0 = assemble((p_i1 - p_0) ** 2 * dx)
        #incr_divu_norm_0 = assemble(div(u_i1 - u_0) ** 2 * dx)

        print("Components of the improved iterative majorants:")
        print("-----------------------------------------------")

        #maj_base_aux = lmbda / 2 * (maj_uhdiv_it[last_it] + maj_uhdiv_it[aux_it]) \
        #               + L / 4 * (maj_phl2_it[last_it] + maj_phl2_it[aux_it])

        # could be a problem?
        incr_sigma_norm_aux = assemble((eta_i1 - eta_aux) ** 2 * dx)
        maj_base_aux = incr_sigma_norm_aux \
                         + lmbda / 2 * (maj_uhdiv_it[last_it] + maj_uhdiv_it[aux_it]) \
                         + L / 4 * (maj_phl2_it[last_it] + maj_phl2_it[aux_it])
        q_aux = q**(test_params["iter_num"] - aux_it)

        print("\| eta_i1 - eta_aux \|  = %.4e" % (incr_sigma_norm_aux))
        #print("maj_uhdiv_it[last_it] + maj_uhdiv_it[aux_it]  = %.4e" % (maj_uhdiv_it[last_it] + maj_uhdiv_it[aux_it]))
        #print("maj_phl2_it[last_it] + maj_phl2_it[aux_it]  = %.4e" % (maj_phl2_it[last_it] + maj_phl2_it[aux_it]))

        print("pow = ", (test_params["iter_num"] - aux_it))
        print("pow = ", (test_params["iter_num"] - aux_it))
        print("q = ", q)
        print("q_aux = q^pow = q^(total-aux) = ", q ** (test_params["iter_num"] - aux_it))
        print("3 * q_aux / (1 - q_aux)**2 * (C_FD ** 2 * beta / min_eig_kappa + 1) = %.4e, " %
              (3 * q_aux / (1 - q_aux) ** 2 * (C_FD ** 2 * funcs["beta"] / funcs["min_eig_kappa"] + 1)))
        #print("where")
        #print("(C_FD ** 2 * beta / min_eig_kappa + 1) = %.4e" % ((C_FD ** 2 * funcs["beta"] / funcs["min_eig_kappa"] + 1)))
        #print("3 * q_0 / (1 - q_0)**2 = %.4e" % (3 * q_aux / (1 - q_aux) ** 2))
        #print("3 * q_0**2 / (1 - q_0)**2 * (1 + gdim * lmbda / (2 * mu)) = %.4e" %(3 * q_aux ** 2 / (1 - q_aux) ** 2 * (1 + gdim * lmbda / (2 * mu))))
        print("maj(i, q_aux) = %.4e\n " % maj_base_aux)

        maj_pi_it_aux = maj_base_aux * 3 * q_aux / (1 - q_aux) ** 2 * (C_FD ** 2 * funcs["beta"] / funcs["min_eig_kappa"] + 1)
        maj_ui_it_aux = maj_base_aux * 3 * q_aux ** 2 / (1 - q_aux) ** 2 * (1 + gdim * lmbda / (2 * mu))

        if test_params["full_documentation"]:
            print("based on q_aux: M^{2,i}(p) = %.4e" % (maj_pi_it_aux / norm_p))
            print("based on q_aux: M^{2,i}(u) = %.4e\n" % (maj_ui_it_aux / norm_u))

        # Scale error and majorant of displacement with tau
        eu_enrg_it[last_it] = self.test_params["tau"] * eu_enrg_it[last_it]
        maj_ui_it[last_it] = self.test_params["tau"] * maj_ui_it[last_it]
        maj_tilde_ui_it[last_it] = self.test_params["tau"] * maj_tilde_ui_it[last_it]
        maj_uh_it[last_it] = self.test_params["tau"] * maj_uh_it[last_it]
        maj_ui_it_aux = self.test_params["tau"] * maj_ui_it_aux

        # Get the value from the last iteration for all the discretisation norms
        ep_enrg, ep_l2, eu_enrg, eu_div, \
        maj_ph, maj_phl2, maj_uh, maj_uhdiv, \
        maj_pi, maj_ui, \
        e, e_p, e_u, \
        maj_p, maj_u, maj, \
        norm_p_accum, norm_u_accum = \
            self.get_last_iter(n, last_it,
                               ep_enrg, ep_l2, eu_enrg, eu_div, maj_ph, maj_phl2, maj_uh, maj_uhdiv, maj_pi, maj_ui,
                               ep_enrg_it, ep_l2_it, eu_enrg_it, eu_div_it,
                               maj_ph_it, maj_phl2_it, maj_uh_it, maj_uhdiv_it,
                               maj_pi_it, maj_ui_it,
                               maj_tilde_pi_it, maj_tilde_ui_it,
                               maj_pi_it_aux, maj_ui_it_aux,
                               e, e_p, e_u, maj_p, maj_u, maj,
                               norm_p, norm_u, norm_p_accum, norm_u_accum)
        '''
        # Get the value from the last iteration for all the discretisation norms
        ep_enrg, ep_l2, eu_enrg, eu_div, maj_ph, maj_phl2, maj_uh, maj_uhdiv, maj_pi, maj_ui, norm_p_accum, norm_u_accum = \
            self.get_last_iter(n, last_it,
                               ep_enrg, ep_l2, eu_enrg, eu_div, maj_ph, maj_phl2, maj_uh, maj_uhdiv, maj_pi, maj_ui,
                               ep_enrg_it, ep_l2_it, eu_enrg_it, eu_div_it,
                               maj_ph_it, maj_phl2_it, maj_uh_it, maj_uhdiv_it,
                               maj_pi_it, maj_ui_it,
                               norm_p, norm_u, norm_p_accum, norm_u_accum)

        # Correction of accumulate majorants with auxuliary step improvement
        # HOW TO SCALE HERE CORRECTLY ??
        if n == 0:
            maj_pi[n] = maj_pi_it_aux / norm_p
            maj_ui[n] = maj_ui_it_aux / norm_u

        else:

            maj_pi[n] = (maj_pi[n - 1] + maj_pi_it_aux) / norm_p_accum[n]
            maj_ui[n] = (maj_ui[n - 1] + maj_ui_it_aux) / norm_u_accum[n]
        '''
        # Print errors and majorants on all time steps
        if test_params["full_documentation"]:
            '''
            print("increment on the %dth interval " % (n+1))
            print("-------------------------------------------------------------------")
            print("||| e_p |||^2  = %.4e, Mp^2  = %.4e, i_eff(Mp)  = %.2f" % (
                  ep_enrg_it[last_it] / norm_p, 2 * (maj_pi_it[last_it] + maj_ph_it[last_it]) / norm_p,
                  sqrt((2 * (maj_pi_it[last_it] + maj_ph_it[last_it])) / ep_enrg_it[last_it])))
            print("||| e_u |||^2  = %.4e, Mu^2  = %.4e, i_eff(Mu)  = %.2f" % (
                  eu_enrg_it[last_it] / norm_u, 2 * (maj_ui_it[last_it] + maj_uh_it[last_it]) / norm_u,
                  sqrt((2 * (maj_ui_it[last_it] + maj_uh_it[last_it])) / eu_enrg_it[last_it])))
            print("-------------------------------------------------------------------")
            print("[(e_u, e_p)]^2 = %.4e, Mit^2 = %.4e, i_eff(Mit) = %.2f"
                  % (ep_enrg_it[last_it] + eu_enrg_it[last_it] / max([norm_p, norm_u]),
                     2 * (maj_pi_it[last_it] + maj_ph_it[last_it] + maj_ui_it[last_it] + maj_uh_it[last_it]) / max([norm_p, norm_u]),
                     sqrt(2 * (maj_pi_it[last_it] + maj_ph_it[last_it] + maj_ui_it[last_it] + maj_uh_it[last_it]) / (ep_enrg_it[last_it] + eu_enrg_it[last_it]))))
            print("-------------------------------------------------------------------")

            print("improved increment on the %dth interval " % (n + 1))
            print("-------------------------------------------------------------------")
            print("||| e_p |||^2  = %.4e, Mp^2  = %.4e, i_eff(Mp)  = %.2f" % (
                ep_enrg_it[last_it] / norm_p,
                2 * (maj_pi_it_aux + maj_ph_it[last_it]) / norm_p,
                sqrt((2 * (maj_pi_it_aux + maj_ph_it[last_it])) / ep_enrg_it[last_it])))
            print("||| e_u |||^2  = %.4e, Mu^2  = %.4e, i_eff(Mu)  = %.2f" % (
                eu_enrg_it[last_it] / norm_u,
                2 * (maj_ui_it_aux + maj_uh_it[last_it]) / norm_u,
                sqrt((2 * (maj_ui_it_aux + maj_uh_it[last_it])) / eu_enrg_it[last_it])))
            print("-------------------------------------------------------------------")
            print("[(e_u, e_p)]^2 = %.4e, Mit^2 = %.4e, i_eff(Mit) = %.2f"
                      % ((ep_enrg_it[last_it] + eu_enrg_it[last_it]) / max([norm_p, norm_u]),
                         2 * (maj_pi_it_aux + maj_ph_it[last_it] + maj_ui_it_aux + maj_uh_it[last_it]) / max([norm_p, norm_u]),
                     sqrt(2 * (maj_pi_it_aux + maj_ph_it[last_it] + maj_ui_it_aux + maj_uh_it[last_it]) / (
                             ep_enrg_it[last_it] + eu_enrg_it[last_it]))))
            print("-------------------------------------------------------------------")
            '''

            print("increment on the %dth interval " % (n+1))
            print("-------------------------------------------------------------------")
            print("||| e_p |||^2  = %.4e, Mp^2  = %.4e, i_eff(Mp)  = %.2f" % (
                  ep_enrg_it[last_it] / norm_p_accum[n],
                  2 * (maj_pi_it[last_it] + maj_ph_it[last_it]) / norm_p_accum[n],
                  sqrt((2 * (maj_pi_it[last_it] + maj_ph_it[last_it])) / ep_enrg_it[last_it])))
            print("||| e_u |||^2  = %.4e, Mu^2  = %.4e, i_eff(Mu)  = %.2f" % (
                  eu_enrg_it[last_it] / norm_u_accum[n],
                  2 * (maj_ui_it[last_it] + maj_uh_it[last_it]) / norm_u_accum[n],
                  sqrt((2 * (maj_ui_it[last_it] + maj_uh_it[last_it])) / eu_enrg_it[last_it])))
            print("-------------------------------------------------------------------")
            print("[(e_u, e_p)]^2 = %.4e, Mit^2 = %.4e, i_eff(Mit) = %.2f"
                  % (ep_enrg_it[last_it] + eu_enrg_it[last_it] / max([norm_p_accum[n], norm_u_accum[n]]),
                     2 * (maj_pi_it[last_it] + maj_ph_it[last_it] + maj_ui_it[last_it] + maj_uh_it[last_it]) / max([norm_p_accum[n], norm_u_accum[n]]),
                     sqrt(2 * (maj_pi_it[last_it] + maj_ph_it[last_it] + maj_ui_it[last_it] + maj_uh_it[last_it]) / (ep_enrg_it[last_it] + eu_enrg_it[last_it]))))
            print("-------------------------------------------------------------------")

            print("improved increment on the %dth interval " % (n + 1))
            print("-------------------------------------------------------------------")
            print("||| e_p |||^2  = %.4e, Mp^2  = %.4e, i_eff(Mp)  = %.2f" % (
                ep_enrg_it[last_it] / norm_p_accum[n],
                2 * (maj_pi_it_aux + maj_ph_it[last_it]) / norm_p_accum[n],
                sqrt((2 * (maj_pi_it_aux + maj_ph_it[last_it])) / ep_enrg_it[last_it])))
            print("||| e_u |||^2  = %.4e, Mu^2  = %.4e, i_eff(Mu)  = %.2f" % (
                eu_enrg_it[last_it] / norm_u_accum[n],
                2 * (maj_ui_it_aux + maj_uh_it[last_it]) / norm_u_accum[n],
                sqrt((2 * (maj_ui_it_aux + maj_uh_it[last_it])) / eu_enrg_it[last_it])))
            print("-------------------------------------------------------------------")
            print("[(e_u, e_p)]^2 = %.4e, Mit^2 = %.4e, i_eff(Mit) = %.2f"
                      % ((ep_enrg_it[last_it] + eu_enrg_it[last_it]) / max([norm_p_accum[n], norm_u_accum[n]]),
                         2 * (maj_pi_it_aux + maj_ph_it[last_it] + maj_ui_it_aux + maj_uh_it[last_it]) / max([norm_p_accum[n], norm_u_accum[n]]),
                     sqrt(2 * (maj_pi_it_aux + maj_ph_it[last_it] + maj_ui_it_aux + maj_uh_it[last_it]) / (
                             ep_enrg_it[last_it] + eu_enrg_it[last_it]))))
            print("-------------------------------------------------------------------")

            print("improved increment (with tilde majorant) on the %dth interval" % (n + 1))
            print("-------------------------------------------------------------------")
            print("||| e_p |||^2  = %.4e, Mp^2  = %.4e, i_eff(Mp)  = %.2f" % (
                ep_enrg_it[last_it] / norm_p_accum[n],
                2 * (maj_tilde_pi_it[last_it] + maj_ph_it[last_it]) / norm_p_accum[n],
                sqrt((2 * (maj_tilde_pi_it[last_it] + maj_ph_it[last_it])) / ep_enrg_it[last_it])))
            print("||| e_u |||^2  = %.4e, Mu^2  = %.4e, i_eff(Mu)  = %.2f" % (
                eu_enrg_it[last_it] / norm_u_accum[n],
                2 * (maj_tilde_ui_it[last_it] + maj_uh_it[last_it]) / norm_u_accum[n],
                sqrt((2 * (maj_tilde_ui_it[last_it]  + maj_uh_it[last_it])) / eu_enrg_it[last_it])))
            print("-------------------------------------------------------------------")
            print("[(e_u, e_p)]^2 = %.4e, Mit^2 = %.4e, i_eff(Mit) = %.2f"
                  % ((ep_enrg_it[last_it] + eu_enrg_it[last_it]) / max([norm_p_accum[n], norm_u_accum[n]]),
                     2 * (maj_tilde_pi_it[last_it] + maj_ph_it[last_it] + maj_tilde_ui_it[last_it]  + maj_uh_it[last_it]) / max(
                         [norm_p_accum[n], norm_u_accum[n]]),
                     sqrt(2 * (maj_tilde_pi_it[last_it] + maj_ph_it[last_it] + maj_tilde_ui_it[last_it]  + maj_uh_it[last_it]) / (
                             ep_enrg_it[last_it] + eu_enrg_it[last_it]))))
            print("-------------------------------------------------------------------\n\n")

            print("-------------------------------------------------------------------")
            print("accumulated result of intervals %d-%dth" % (0, n + 1))
            print("-------------------------------------------------------------------")
            print("||| e_p |||^2  = %.4e, Mp^2 = %.4e, i_eff(Mp) = %.2f" % (
                e_p[n] / norm_p_accum[n],
                maj_p[n] / norm_p_accum[n],
                sqrt(maj_p[n] / e_p[n])))
            print("||| e_u |||^2  = %.4e, Mu^2 = %.4e, i_eff(Mu) = %.2f" % (
                e_u[n] / norm_u_accum[n],
                maj_u[n] / norm_u_accum[n],
                sqrt(maj_u[n] / e_u[n])))
            print("-------------------------------------------------------------------")
            print("[(e_u, e_p)]^2 = %.4e, M^2  = %.4e, i_eff(M) = %.2f" %
                  (e[n] / max([norm_p_accum[n], norm_u_accum[n]]),
                   maj[n] / max([norm_p_accum[n], norm_u_accum[n]]),
                   sqrt(maj[n] / e[n])))
            print("-------------------------------------------------------------------")
        '''
        # Accumulate errors and majorants on all time steps
        e_p[n] = ep_enrg[n]
        e_u[n] = eu_enrg[n]
        e[n] = ep_enrg[n] + eu_enrg[n]

        maj_p[n] = 2 * (maj_pi[n] + maj_ph[n])
        maj_u[n] = 2 * (maj_ui[n] + maj_uh[n])
        maj[n]   = (maj_p[n] + maj_u[n])

        if test_params["full_documentation"]:
            print("accumulated result of intervals %d-%dth" % (0, n+1))
            print("-------------------------------------------------------------------")
            print("||| e_p |||^2  = %.4e, Mp^2 = %.4e, i_eff(Mp) = %.2f" % (
                e_p[n], maj_p[n], sqrt(maj_p[n] / e_p[n])))
            print("||| e_u |||^2  = %.4e, Mu^2 = %.4e, i_eff(Mu) = %.2f" % (
                e_u[n], maj_u[n], sqrt(maj_u[n] / e_u[n])))
            print("-------------------------------------------------------------------")
            print("[(e_u, e_p)]^2 = %.4e, M^2  = %.4e, i_eff(M) = %.2f" % (e[n], maj[n], sqrt(maj[n] / e[n])))
            print("-------------------------------------------------------------------")
        #
        '''
        #input("Press")
        return u_i1, p_i1, \
               ep_enrg, ep_l2, eu_enrg, eu_div, \
               e_p, e_u, e, \
               maj_ph, maj_phl2, maj_uh, maj_uhdiv, maj_pi, maj_ui, \
               maj_p, maj_u, maj, \
               norm_p_accum, norm_u_accum

    # Function with full-implicit procedure
    def full_implicit_coupling(self,
                               n, bc_p, bc_u,
                               alpha, beta, g_tilde,
                               mu, lmbda,
                               error_control, ep_enrg, ep_l2, eu_enrg, eu_div, e_enrg, maj,
                               func_spaces, funcs):

        ep_enrg_it = error_control.get_error_array(1)
        eu_enrg_it = error_control.get_error_array(1)
        ep_l2_it = error_control.get_error_array(1)
        eu_div_it = error_control.get_error_array(1)
        maj_pd = error_control.get_error_array(self.test_params["time_steps"])
        maj_peq = error_control.get_error_array(self.test_params["time_steps"])
        maj_ud = error_control.get_error_array(self.test_params["time_steps"])
        maj_ueq = error_control.get_error_array(self.test_params["time_steps"])

        (u, p) = TrialFunctions(func_spaces["W"])
        (v, theta) = TestFunctions(func_spaces["W"])
        w_n1 = Function(func_spaces["W"])

        # Collect BC for the mixed method
        bcs = [bc_u, bc_p]

        epsilon = lambda w: sym(grad(w))
        sigma   = lambda w: 2 * mu * epsilon(w) + lmbda * div(w) * Identity(self.domain_params["gdim"])#tr(epsilon(w)) * Identity(self.domain_params["gdim"])

        # Variational forms of the Biot RHS and LHS
        biot_lhs = (inner(beta * p, theta) + inner(funcs["kappa"] * grad(p), grad(theta))
                    + alpha * inner(div(u), theta)
                    + tr(inner(sigma(u), epsilon(v)))
                    - alpha * inner(p, div(v))) * dx
        biot_rhs = (inner(funcs["f"], v) + inner(g_tilde, theta)) * dx

        # Solve the Biot system
        solve(biot_lhs == biot_rhs, w_n1, bcs)
        u_n1, p_n1 = w_n1.split()

        # Calculate the error in the pressure term
        ep_enrg_it, ep_l2_it, var_grad_e_p, var_e_p, norm_p  = error_control.error_pressure(ep_enrg_it, ep_l2_it, 1, p_n1, funcs)
        eu_enrg_it, eu_div_it, var_strain_eu, var_div_eu, norm_u = error_control.error_displacement(eu_enrg_it, eu_div_it, 1, u_n1, funcs, mu, lmbda)

        funcs["react"] = beta
        G = g_tilde - alpha * div(u_n1) - beta * p_n1
        F = funcs["f"] - alpha * grad(p_n1)

        maj_pd, maj_peq, var_md_ep = \
            error_control.majorant_pressure_part(n, maj_pd, maj_peq, p_n1, G,
                                            self.domain_params, self.test_params,
                                            func_spaces, funcs)
        maj_ud, maj_ueq, var_md_eu = \
            error_control.majorant_displacement_part(n, maj_ud, maj_ueq, u_n1, F, mu, lmbda,
                                                 self.domain_params, self.test_params,
                                                 func_spaces, funcs)
        
        if n == 0 :
            ep_enrg[n] = ep_enrg_it[0]
            ep_l2[n]   = ep_l2_it[0]
            eu_enrg[n] = eu_enrg_it[0]
            eu_div[n]  = eu_div_it[0]
            e_enrg[n]  = ep_enrg[n] + eu_enrg[n]
            maj[n]     = maj_pd[n] + maj_peq[n] + maj_ud[n] + maj_ueq[n]
        else:
            ep_enrg[n] = ep_enrg[n-1] + ep_enrg_it[0]
            ep_l2[n]   = ep_l2[n-1] + ep_l2_it[0]
            eu_enrg[n] = eu_enrg[n-1] + eu_enrg_it[0]
            eu_div[n]  = eu_div[n-1] + eu_div_it[0]

            e_enrg[n]  = e_enrg[n-1] + ep_enrg[n] + eu_enrg[n]
            maj[n]     = maj[n-1] + maj_pd[n] + maj_peq[n] + maj_ud[n] + maj_ueq[n]

        return u_n1, p_n1, ep_enrg, ep_l2, eu_enrg, eu_div, e_enrg, maj


# Class of the Dirichlet BC
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


# Class for testing Biot problem
class TestBiot():
    # Init variable of TestBiot class
    def __init__(self, problem_data, domain_params, material_params, test_params):

        self.problem_data = problem_data
        self.domain_params = domain_params
        self.material_params = material_params
        self.test_params = test_params

    def get_initial_pressure_and_displacement(self, ds, n, p, u, v, theta, bc_p, bc_u, p_k, u_k, f, t_N, phi_f, g, alpha, mu, lmbda):

        # Define the variational formulation for pressure
        a_p0 = inner(grad(p), grad(theta)) * dx
        l_p0 = phi_f * inner(Constant((0, -g)), grad(theta)) * dx

        # Assemble the system for the initial pressure
        A_p0, b_p0 = assemble_system(a_p0, l_p0, bc_p)
        solve(A_p0, p_k.vector(), b_p0)

        epsilon = lambda v: 0.5 * (nabla_grad(v) + nabla_grad(v).T)
        sigma = lambda v: 2 * mu * epsilon(v) + lmbda * div(v) * Identity(self.self.domain_params['gdim'])

        # Define the variationalphi_0 formulation for displacement
        a_u0 = tr(inner(sigma(u), epsilon(v))) * dx
        l_u0 = inner(f - alpha * grad(p_k), v) * dx + inner(t_N, v) * ds(4)

        # Assemble the system for initial displacement
        A_u0, b_u0 = assemble_system(a_u0, l_u0, bc_u)
        solve(A_u0, u_k.vector(), b_u0)

        return p_k, u_k

    # Construction of the required functional space
    def functional_spaces(self, mesh, test_params):

        # Order of the space of the exact solutions
        v_exact_degree = test_params['u_approx_order'] + 2
        p_exact_degree = test_params['p_approx_order'] + 2

        # Space for the displacement
        V       = VectorFunctionSpace(mesh, "CG", test_params['u_approx_order'])
        V_EL    = VectorElement("CG", mesh.ufl_cell(), test_params['u_approx_order'])
        V_exact = VectorFunctionSpace(mesh, "CG", v_exact_degree)
        Vex_EL  = VectorElement("CG", mesh.ufl_cell(), v_exact_degree)

        # Space for the pressure
        Q       = FunctionSpace(mesh, "CG", test_params['p_approx_order'])
        Q_EL    = FiniteElement("CG", mesh.ufl_cell(), test_params['p_approx_order'])
        Q_exact = FunctionSpace(mesh, "CG", p_exact_degree)
        Qex_EL  = FiniteElement("CG", mesh.ufl_cell(), p_exact_degree)

        W = FunctionSpace(mesh, V_EL * Q_EL)
        W_exact = FunctionSpace(mesh, Vex_EL * Qex_EL)

        H_div = FunctionSpace(mesh, 'RT', test_params['flux_approx_order'])
        H_Div = TensorFunctionSpace(mesh, 'CG', test_params['stress_approx_order'])

        func_spaces = dict(V=V,
                           Q=Q,
                           W=W,
                           V_exact=V_exact,
                           Q_exact=Q_exact,
                           H_div=H_div,
                           H_Div=H_Div,
                           W_exact=W_exact)
        return func_spaces

    def update_term_dependent_on_time(self, funcs, t):

        if self.test_params['test_num'] == 3:
            funcs['f'].t = t
            funcs['p_D_gamma1'].t = t
            funcs['g'].t = t

        else:
            funcs['u_e'].t = t
            funcs['p_e'].t = t
            funcs['f'].t = t
            funcs['g'].t = t

        return funcs

    def convert_expressions_to_functions(self, func_spaces):
        '''
        :return:
        '''
        if self.test_params['test_num'] == 1:

            # Exact pressure
            p_e = Expression(self.problem_data['p_expr'],
                             F=self.material_params['F'],
                             B=self.material_params['B'],
                             nu_u=self.material_params['nu_u'],
                             a=self.domain_params['l_x'],
                             alpha_1=self.material_params['alpha_n'][0],
                             alpha_2=self.material_params['alpha_n'][1],
                             alpha_3=self.material_params['alpha_n'][2],
                             alpha_4=self.material_params['alpha_n'][3],
                             alpha_5=self.material_params['alpha_n'][4],
                             c_f=self.material_params['c_f'],
                             t=0.0,
                             degree=4,
                             element=func_spaces["Q_exact"].ufl_element())
            # Initial condition for pressure
            p_0 = p_e
            # Dirichlet BC for pressure
            # p_D = Expression(self.problem_data['pD_expr'])
            # Flow source/sink
            g = Expression(self.problem_data['g_expr'], degree=4)
            f = Expression((self.problem_data['f_expr'][0], self.problem_data['f_expr'][1]), degree=4)
            u_e = Expression((self.problem_data['u_expr'][0], self.problem_data['u_expr'][1]),
                             F=self.material_params['F'],
                             nu_u=self.material_params['nu_u'],
                             nu=self.material_params['nu'],
                             a=self.domain_params['l_x'],
                             mu=self.material_params['mu'],
                             alpha_1=self.material_params['alpha_n'][0],
                             alpha_2=self.material_params['alpha_n'][1],
                             alpha_3=self.material_params['alpha_n'][2],
                             alpha_4=self.material_params['alpha_n'][3],
                             alpha_5=self.material_params['alpha_n'][4],
                             c_f=self.material_params['c_f'],
                             t=0.0,
                             degree=4)
            u_e_0 = Expression(self.problem_data['u_expr'][0],
                               F=self.material_params['F'],
                               nu_u=self.material_params['nu_u'],
                               nu=self.material_params['nu'],
                               a=self.domain_params['l_x'],
                               mu=self.material_params['mu'],
                               alpha_1=self.material_params['alpha_n'][0],
                               alpha_2=self.material_params['alpha_n'][1],
                               alpha_3=self.material_params['alpha_n'][2],
                               alpha_4=self.material_params['alpha_n'][3],
                               alpha_5=self.material_params['alpha_n'][4],
                               c_f=self.material_params['c_f'],
                               t=0.0,
                               degree=4)
            u_e_1 = Expression(self.problem_data['u_expr'][1],
                               F=self.material_params['F'],
                               nu_u=self.material_params['nu_u'],
                               nu=self.material_params['nu'],
                               a=self.domain_params['l_x'],
                               mu=self.material_params['mu'],
                               alpha_1=self.material_params['alpha_n'][0],
                               alpha_2=self.material_params['alpha_n'][1],
                               alpha_3=self.material_params['alpha_n'][2],
                               alpha_4=self.material_params['alpha_n'][3],
                               alpha_5=self.material_params['alpha_n'][4],
                               c_f=self.material_params['c_f'],
                               t=0.0, degree=4)
            t_N = Expression((self.problem_data['t_N_expr'][0],
                              self.problem_data['t_N_expr'][1]),
                              F=self.material_params['F'], degree=4)
            u_0 = u_e

            funcs = dict(u_e=u_e,
                         u_e_0=u_e_0,
                         u_e_1=u_e_1,
                         p_e=p_e,
                         u_0=u_0,
                         p_0=p_0,
                         f=f,
                         t_N=t_N,
                         g=g)

        elif self.test_params['test_num'] == 2 or \
                self.test_params['test_num'] == 4 or \
                self.test_params['test_num'] == 101 or \
                self.test_params['test_num'] == 102 or \
                self.test_params['test_num'] == 103 or \
                self.test_params['test_num'] == 104:
            # Exact pressure and displacement
            p_e = Expression(self.problem_data['p_expr'], t=0.0, degree=4)
            u_e = Expression((self.problem_data['u_expr'][0], self.problem_data['u_expr'][1]), t=0, degree=4)

            # Initial condition for pressure and displacement
            p_0 = Expression(self.problem_data['p0_expr'],degree=4)
            u_0 = Expression((self.problem_data['u0_expr'][0], self.problem_data['u0_expr'][1]), degree=4)

            # Flow source/sink
            g = Expression(self.problem_data['g_expr'],
                           beta=self.material_params['beta'],
                           alpha=self.material_params['alpha'],
                           t=0.0, degree=4)

            # Mechanical load
            f = Expression((self.problem_data['f_expr'][0], self.problem_data['f_expr'][1]),
                           mu=self.material_params['mu'],
                           lmbda=self.material_params['lmbda'],
                           alpha=self.material_params['alpha'],
                           t=0.0, degree=4)
            t_N = Expression((self.problem_data['t_N_expr'][0], self.problem_data['t_N_expr'][1]), degree=4)

            funcs = dict(u_e=u_e,
                         p_e=p_e,
                         u_0=u_0,
                         p_0=p_0,
                         f=f,
                         g=g,
                         t_N=t_N)

        elif self.test_params['test_num'] == 3:
            # Dirichlet boundary conditions for pressure
            p_D_gamma1 = Expression(self.problem_data['p_D_expr'][0], t=0.0, degree=4)
            p_D_gamma = Expression(self.problem_data['p_D_expr'][1], degree=4)

            # Initial condition for pressure
            p_0 = Expression(self.problem_data['p0_expr'], degree=4)

            # Neumann boundary conditions for pressure
            t_N = Constant((self.problem_data['t_N_expr'][0], self.problem_data['t_N_expr'][1]), degree=4)
            u_0 = Expression((self.problem_data['u0_expr'][0], self.problem_data['u0_expr'][1]), degree=4)

            g = Expression(self.problem_data['g_expr'], degree=3)
            f = Expression(self.problem_data['f_expr'][0], self.problem_data['f_expr'][1], degree=3)

            funcs = dict(p_D_gamma1=p_D_gamma1,
                         p_D_gamma=p_D_gamma,
                         u_0=u_0,
                         p_0=p_0,
                         f=f,
                         g=g,
                         t_N=t_N)
        return funcs

    def test_biot(self, test_params, project_path):

        # Check if exist and create folder for the results
        resultfolder_name = postprocess.construct_result_folder_name(self.domain_params['gdim'],
                                                                     self.test_params['test_num'], test_params)
        postprocess.create_results_folder(test_params, resultfolder_name)

        # Generate mesh and boundary parts
        if self.domain_params['gdim'] == 2:
            mesh, facet_function, h = problem.geometry_2d(self.domain_params, self.test_params['test_num'],
                                                          test_params["mesh_resolution"])
        elif self.domain_params['gdim'] == 3:
            mesh, facet_function, h = problem.geometry_3d(self.domain_params, self.test_params['test_num'])
        # Get the mesh from facet_function
        # mesh = facet_function.mesh()

        '''
        if self.domain_params['gdim'] == 2:
            postprocess.plot_mesh(mesh, project_path + resultfolder_name + 'initial-mesh')
        elif self.domain_params['gdim'] == 3:
            postprocess.plot_mesh_3d(mesh, project_path + resultfolder_name + 'initial-mesh')
        '''
        # Construct reuired functional spaces
        func_spaces = self.functional_spaces(mesh, test_params)

        # Construct required functions
        funcs = self.convert_expressions_to_functions(func_spaces)

        # Set Dirichlet BC for the problem,
        # bc_p, bc_u = self.dirichlet_boundary_conditions(func_spaces, funcs, facet_function, material_params)

        postprocess.output_problem_characteristics(test_num, mesh, problem_data, domain_params, material_params,
                                                   test_params, func_spaces)

        # Set start time
        t0_problem = process_time()

        # Solve Biot system
        error_p, errpr_pl2, error_u, error_udiv, e_p, e_u, e_enrg, maj_p, maj_u, maj, accum_norm_p, accum_norm_u \
            = self.solve_biot(func_spaces, funcs, facet_function, mesh,
                              test_params, project_path, resultfolder_name)
        '''
        ep_enrg[test_params["time_steps"] - 1], ep_l2[test_params["time_steps"] - 1], \
        eu_enrg[test_params["time_steps"] - 1], eu_div[test_params["time_steps"] - 1], \
        maj_p[test_params["time_steps"] - 1], maj_u[test_params["time_steps"] - 1], \
        maj[test_params["time_steps"] - 1]
        '''
        # Set the final time
        tfinal_problem = process_time()
        print('total values:')
        print('\t\th\t||| e_p |||^2\t\t\tMp^2\t||| e_u |||^2\t\t\tMu^2\t[(e_u, e_p)]^2\t M^2\t i_eff')
        print('%4.4e &\t %12.4e &\t %12.4e &\t %12.4e &\t %12.4e &\t %12.4e &\t%12.4e &\t%5.2f\n'
              % (h,
                 e_p[-1] / accum_norm_p[-1],
                 maj_p[-1] / accum_norm_p[-1],
                 e_u[-1] / accum_norm_u[-1],
                 maj_u[-1] / accum_norm_u[-1],
                 e_enrg[-1] / max([accum_norm_p[-1], accum_norm_u[-1]]),
                 maj[-1] / max([accum_norm_p[-1], accum_norm_u[-1]]),
                 sqrt(maj[test_params["time_steps"] - 1] / e_enrg[test_params["time_steps"] - 1])))
        '''
        print('%4.4e &\t %12.4e &\t %12.4e &\t %12.4e &\t %12.4e &\t %12.4e &\t%12.4e &\t%5.2f\n'
              % (h,
                 e_p[test_params["time_steps"] - 1],
                 maj_p[test_params["time_steps"] - 1],
                 e_u[test_params["time_steps"] - 1],
                 maj_u[test_params["time_steps"] - 1],
                 e_enrg[test_params["time_steps"] - 1],
                 maj[test_params["time_steps"] - 1],
                 sqrt(maj[test_params["time_steps"] - 1] / e_enrg[test_params["time_steps"] - 1])))
    
        '''
        '''
        print('\th\t||| e_p |||^2\t|| e_p ||^2\t||| e_u |||^2\t|| div(u) ||^2\t\t\tMp^2\t\t\tMu^2\ti_eff')
        print('%4.4e\t %12.4e\t %10.4e\t %12.4e\t %13.4e %13.4e\t %11.4e\t%5.2f\n'
              % (h,
                 error_p[test_params["time_steps"] - 1], errpr_pl2[test_params["time_steps"] - 1],
                 error_u[test_params["time_steps"] - 1], error_udiv[test_params["time_steps"] - 1],
                 maj_p[test_params["time_steps"] - 1], maj_u[test_params["time_steps"] - 1],
                 sqrt(maj[test_params["time_steps"] - 1] / e_enrg[["time_steps"] - 1])))
        '''
        '''
        print('total values')
        print('\t\th\t||| e_p |||^2\t|| e_p ||^2\t||| e_u |||^2\t|| div(u) ||^2\t\t\tMp^2\t\t\tMu^2\ti_eff')
        print('%4.4e\t %12.4e\t %10.4e\t %12.4e\t %13.4e %13.4e\t %11.4e\t%5.2f\n'
              % (h, sum(error_p), sum(errpr_pl2), sum(error_u), sum(error_udiv),
                 sum(maj_p), sum(maj_u), sqrt(sum(maj) / (sum(error_p) + sum(error_u)))))

        # print("total time = %.2f secs" % (tfinal_problem - t0_problem))
        print("total time = %.2f mins" % (float(tfinal_problem - t0_problem) / 60.0))
        '''
        # interactive()

    def fetch_paramenters_with_scaling(self, tau):

        # Flow-equation realted parameters
        alpha = float(self.material_params["alpha"])
        beta = float(self.material_params["beta"])
        mu_f = float(self.material_params["mu_f"])
        # Getting kappa, inv_kappa, min_eig_kappa with scaling
        min_eig_kappa = tau / mu_f * float(self.material_params["min_eig_kappa"])
        # Convert matrix to UFL form
        kappa = tau / mu_f * as_matrix(self.material_params["kappa"])

        # Mechanics-equation realted parameters
        mu = float(self.material_params["mu"])
        lmbda = float(self.material_params["lmbda"])

        return alpha, beta, kappa, min_eig_kappa, lmbda, mu

    def solve_biot(self, func_spaces, funcs, facet_function, mesh, test_params, project_path, resultsfolder_name):

        # Define the value of the time-step
        tau = float(self.problem_data['t_T'] / test_params['time_steps'])
        self.test_params['tau'] = tau

        # Initialize times
        t_n = 0  # t_n
        t_n1 = tau  # t_{n+1}

        # Initialial time counter
        n = 0

        # Fetch required parameters with scaling by tau
        alpha, beta, kappa, min_eig_kappa, lmbda, mu \
            = self.fetch_paramenters_with_scaling(tau)

        funcs["beta"] = beta
        funcs["kappa"] = kappa
        funcs["min_eig_kappa"] = min_eig_kappa

        # Geometry parameters
        ds = Measure("ds", subdomain_data=facet_function)
        normal = FacetNormal(mesh)
        gdim = self.domain_params['gdim']

        # Define function p_{n-1} and u_{n-1}, which is p_0 and u_0 on the initial step
        u_n = interpolate(funcs['u_0'], func_spaces['V'])
        p_n = interpolate(funcs['p_0'], func_spaces['Q'])

        # Define the error control class
        error_control = ErrorControl(mesh,
                                     func_spaces['V_exact'], func_spaces['Q_exact'],
                                     self.domain_params["gdim"])

        norm_p_accum = error_control.get_error_array(test_params["time_steps"])
        norm_u_accum = error_control.get_error_array(test_params["time_steps"])

        ep_enrg = error_control.get_error_array(test_params["time_steps"])
        eu_enrg = error_control.get_error_array(test_params["time_steps"])
        ep_l2 = error_control.get_error_array(test_params["time_steps"])
        eu_div = error_control.get_error_array(test_params["time_steps"])
        ep = error_control.get_error_array(test_params["time_steps"])
        eu = error_control.get_error_array(test_params["time_steps"])
        e = error_control.get_error_array(test_params["time_steps"])

        maj_ph = error_control.get_error_array(test_params["time_steps"])
        maj_phl2 = error_control.get_error_array(test_params["time_steps"])
        maj_uh = error_control.get_error_array(test_params["time_steps"])
        maj_uhdiv = error_control.get_error_array(test_params["time_steps"])
        maj_pi = error_control.get_error_array(test_params["time_steps"])
        maj_ui = error_control.get_error_array(test_params["time_steps"])
        maj_p = error_control.get_error_array(test_params["time_steps"])
        maj_u = error_control.get_error_array(test_params["time_steps"])
        maj = error_control.get_error_array(test_params["time_steps"])

        # Define the boit solver class
        biot_solver = BiotSolvers(mesh, self.material_params, self.domain_params, self.test_params)

        # Loop over time-steps
        while t_n1 <= self.problem_data['t_T'] + DOLFIN_EPS:

            print('%----------------------------------------------------')
            print(' Time step (t_%i, t_%i) = (%.3f, %.3f), tau = %.2e:' % (n, n+1, t_n, t_n1, tau))
            print('%----------------------------------------------------')

            # Update the fource term
            funcs = self.update_term_dependent_on_time(funcs, t_n1)

            # Define the source function on the n-th time step
            g_tilde = tau * funcs["g"] + beta * p_n + alpha * div(u_n)#tr(epsilon(u_n))#

            # Update Dirichlet BC if the given data is dependent on time
            bc_p, bc_u = self.dirichlet_boundary_conditions(func_spaces, funcs, facet_function, material_params)

            if test_params['coupling_approach'] == 'fully-implicit':

                u_n1, p_n1, ep_enrg, ep_l2, eu_enrg, eu_div, e, maj \
                    = biot_solver.full_implicit_coupling(n, bc_p, bc_u, alpha, beta, g_tilde,
                                                         mu, lmbda,
                                                         error_control, ep_enrg, ep_l2, eu_enrg, eu_div, e, maj,
                                                         func_spaces, funcs)
            else:

                # Contractive scheme parameters
                L = alpha ** 2 / (2 * lmbda)
                L_opt = alpha ** 2 / (2 * (lmbda + 2 * mu / gdim))
                q = L / (L + beta)
                q_opt = L_opt / (L_opt + beta)

                funcs["react"] = L_opt + beta
                postprocess.print_biot_iterative_coupling_parameters(beta, L, L_opt, q, q_opt)

                # Define trial fucntions for mechanics and flow equations
                u = TrialFunction(func_spaces['V'])
                p = TrialFunction(func_spaces['Q'])

                # Define test functions for mechanics and flow equations
                v = TestFunction(func_spaces['V'])
                theta = TestFunction(func_spaces['Q'])

                # Define function p_{n, k},p_{n, k+1} and u_{n, k}, u_{n, k+1}
                u_k1 = Function(func_spaces['V'])
                p_k1 = Function(func_spaces['Q'])

                # Set initial values of pressure and displacement in iteration procedure
                p_k = interpolate(Expression("0", degree=3), func_spaces['Q'])
                u_k = interpolate(Expression(("0", "0"), degree=3), func_spaces['V'])
                # p_k = interpolate(funcs["p_e"], func_spaces['Q'])
                # u_k = interpolate(funcs["u_e"], func_spaces['V'])
                # u_k = u_n
                # p_k = p_n
                """
                p_k, u_k = self.get_initial_pressure_and_displacement(ds, normal,
                                                                      p, u, v, theta, bc_p, bc_u,
                                                                      p_k, u_k,
                                                                      funcs['f'], funcs['t_N'],
                                                                      0.2, 9.8, alpha, mu, lmbda)
                """
                u_n1, p_n1, \
                ep_enrg, ep_l2, eu_enrg, eu_div, ep, eu, e, \
                maj_ph, maj_phl2, maj_uh, maj_uhdiv, \
                maj_pi, maj_ui, \
                maj_p, maj_u, maj , \
                norm_p_accum, norm_u_accum = biot_solver.iterative_coupling(n, u, p, v, theta, bc_p, bc_u,
                                                     error_control, ep_enrg, ep_l2, eu_enrg, eu_div, ep, eu, e,
                                                     maj_ph, maj_phl2, maj_uh, maj_uhdiv, maj_pi, maj_ui, maj_p, maj_u,
                                                     maj,
                                                     alpha, mu, lmbda, L_opt, q_opt, ds,
                                                     u_k, u_k1, p_k, p_k1,
                                                     g_tilde,
                                                     funcs, func_spaces, test_params,
                                                     norm_p_accum, norm_u_accum)
            # Update time layer
            t_n1 += tau
            t_n += tau
            n += 1

            if t_n < self.problem_data['t_T'] + DOLFIN_EPS and t_n1 >= self.problem_data['t_T']:
                # Plot pressure
                postprocess.plot_function(p_n1, mesh, gdim,
                                          project_path + resultsfolder_name + 'pressure-n-%d' % (n+1))
                # Plot displacement
                postprocess.plot_vector_function(u_n1, mesh, gdim,
                                                 project_path + resultsfolder_name + 'displacement-n-%d' % (n+1))

                '''
                postprocess.plot_function(interpolate(funcs["p_e"], func_spaces["Q_exact"]),
                                          mesh, gdim,
                                          project_path + resultsfolder_name + 'pressure-exact-n-%d' % (n+1))
                postprocess.plot_vector_function(interpolate(funcs["u_e"], func_spaces["V_exact"]),
                                                 mesh, gdim,
                                                 project_path + resultsfolder_name + 'displacement-exact-n-%d' % (
                                                         n+1))
                '''
            # Update functions
            if self.test_params["coupling_approach"] == "fully-implicit":
                #u_n = interpolate(funcs['u_e'], func_spaces['V'])
                #p_n = interpolate(funcs['p_e'], func_spaces['Q'])
                assign(u_n.sub(0), u_n1.sub(0))
                assign(u_n.sub(1), u_n1.sub(1))
                assign(p_n, p_n1)

            elif self.test_params["coupling_approach"] == "iterative":
                #u_n = interpolate(funcs['u_e'], func_spaces['V'])
                #p_n = interpolate(funcs['p_e'], func_spaces['Q'])
                assign(u_n, u_n1)   # u_{n} <- u_{n+1}
                assign(p_n, p_n1)   # p_{n} <- p_{n+1}
                #u_n.assign(u_n1)
                #p_n.assign(p_n1)
            '''
            postprocess.output_errors_and_estimates_wrt_times(ep_enrg_it, ep_l2_it, eu_enrg_it, eu_div_it,
                                                       maj_ph_it, maj_phl2_it, maj_uh_it, maj_uhdiv_it,
                                                       maj_pi_it, maj_ui_it,
                                                       test_params["iter_num"])
            '''
        return ep_enrg, ep_l2, eu_enrg, eu_div, ep, eu, e, \
               maj_p, maj_u, maj, norm_p_accum, norm_u_accum

    def dirichlet_boundary_conditions(self, func_spaces, funcs, facet_function, material_params):

        # left: 1
        # right: 2
        # bottom: 3
        # top: 4

        if self.test_params['test_num'] == 1:
            # funcs['u_D'], funcs['p_D'],
            # Set Dirichlet BC for pressure on the right(2) boundary of the domain
            bcs_p = [DirichletBC(func_spaces['Q'], Constant(0.0), facet_function, 2)]

            # Set Dirichlet BC for displacement on the left(1) and bottom(3) boundary of the domain
            bcs_u = [DirichletBC(func_spaces['V'].sub(0), Constant(0.0), facet_function, 1),  # u_x = 0
                     DirichletBC(func_spaces['V'].sub(1), Constant(0.0), facet_function, 3),  # u_y = 0
                     DirichletBC(func_spaces['V'].sub(1), funcs["u_e_1"], facet_function, 4)]  # u_y = U_2
            """
            bcs_u = [DirichletBC(func_spaces['V'].sub(0), Constant(0.0), facet_function, 1),
                     DirichletBC(func_spaces['V'].sub(1), Constant(0.0), facet_function, 3),
                     DirichletBC(func_spaces['V'].sub(1), Constant(2*material_params["F"]), facet_function, 4)]
            """
        elif self.test_params['test_num'] == 2 or \
                self.test_params['test_num'] == 4 or \
                self.test_params['test_num'] == 101 or \
                self.test_params['test_num'] == 102 or \
                self.test_params['test_num'] == 103 or \
                self.test_params['test_num'] == 104:
            # Set Dirichlet BC for pressure and displacement on the whole boundary
            if self.test_params["coupling_approach"] == "fully-implicit":
                bcs_p = DirichletBC(func_spaces["W"].sub(1), funcs['p_e'], DirichletBoundary())
                bcs_u = DirichletBC(func_spaces["W"].sub(0), funcs['u_e'], DirichletBoundary())
            elif self.test_params["coupling_approach"] == "iterative":
                bcs_p = DirichletBC(func_spaces['Q'], funcs['p_e'], DirichletBoundary())
                bcs_u = DirichletBC(func_spaces['V'], funcs['u_e'], DirichletBoundary())

        elif self.test_params['test_num'] == 3:
            bcs_p = [DirichletBC(func_spaces['Q'], funcs['p_e'], facet_function, 1),
                     DirichletBC(func_spaces['Q'], funcs['p_e'], facet_function, 0)]
            bcs_u = []

        return bcs_p, bcs_u


if __name__ == '__main__':
    # Defnie the file path and project path
    project_path = postprocess.get_project_path()

    # Types of domain
    rectangular_domain_tag = "rectangular-domain"
    rectangular_domain_with_obsticle_tag = "rectangular-with-obstacle"

    # Decoupling
    fully_implicit_coupling = 'fully-implicit'
    iterative_coupling = 'iterative'

    # Error format
    relative_error = "relative_error"
    absolute_error = "absolute_error"

    # Refinement
    uniform = "uniform"
    adaptive = "adaptive"

    # Pressure recovery method
    CG_method = 'CG-method'
    mixed_method = 'mixed-method'

    # Pressure recovery method
    console_output = 'console'
    file_output = 'file'

    # list of tests
    tests = {1: tests.mandels_problem_2d_t,
             2: tests.simple_example_2d_t,
             3: tests.kolesov_example_2d_t,
             4: tests.simple_example_2_2d_t, # example 1 (in the CAM paper)
             5: tests.simple_example_3d_t,
             101: tests.simple_example_2_2d_t_EGPa,
             102: tests.simple_example_2_2d_t_EGPa_K100,
             103: tests.simple_example_2d_t_bothetall_paper_parameters,
             104: tests.simple_example_2d_t_lubrication_paper_parameters}
    # Set the number of the test and call for the problem data
    #test_num = 102
    #test_num = 2
    test_num = 104
    #test_num = 102

    #resolutions = [16]
    #resolutions = [64]
    #resolutions = [32]
    #resolutions = [64]
    #resolutions = [4]
    #resolutions = [8]
    #resolutions = [8, 16, 32, 64]
    #resolutions = [8, 16, 32, 64]

    resolutions = [4, 8, 16, 32, 64]
    #resolutions = [4, 8, 16, 32, 64]
    #resolutions = [16, 32, 64]
    #resolutions = [64]
    #resolutions = [16]

    for i in range(0, len(resolutions)):
        # Init problem parameters
        test_params = dict(u_approx_order=1,  # Functional spaces parameters for displacement
                           p_approx_order=1,  # Functional spaces parameters for pressure
                           flux_approx_order=2,
                           stress_approx_order=2,
                           iter_accuracy=1e-4,  # Required accuracy at each interation cycle
                           time_steps=100,  # Number of time steps on the interval [0, t_T]
                           mesh_resolution=resolutions[i],  # Lever of refinement of initial mesh [4, 8, 16, 32, 64, 128]
                           iter_num=5,
                           pow=3,
                           coupling_approach=iterative_coupling,
                           pressure_recovery_method=CG_method,
                           full_documentation=True,
                           error_format=relative_error,
                           error_estimates=True,
                           majorant_optimisation=True,
                           majorant_optimization_iterations_number=0,
                           test_num=test_num,
                           output=file_output)

        problem_data, domain_params, material_params = tests[test_num]()

        test = TestBiot(problem_data, domain_params, material_params, test_params)
        test.test_biot(test_params, project_path)
        
