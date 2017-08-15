__author__ = "Svetlana Matculevich <svetlana.v.matculevich@jyu.fi>"
__date__ = "2016-04-20"
__copyright__ = "Copyright (C) 2016 Svetlana Matculevich"
__license__ = ""

from dolfin import *
from dolfin.cpp.io import File, interactive
from dolfin.cpp.mesh import SubDomain
import time
import numpy as np

# Customer libraries
import tests, problem, postprocess
from problem import Grad, Div, strain_tensor

parameters["plotting_backend"] = "matplotlib"

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
        return np.zeros(num)

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

    def error_pressure(self, error_array, i, p, pe, mu_f, kappa, tau, beta):
        """
        :param p: approximate pressure
        :param pe: exact presure
        :param mu_f: fluid ...
        :param kappa: permiability tensor
        :param tau: time step
        :param beta: (1/M + phi_0 * c_f) fluid characteristic
        :param Ve: functional space of exact solution
        :param mesh: mesh
        :param dim: dimension of the problem

        :return: error-norm between p and ue
        """
        # Interpolate exact and approximate solution to the functional space of exact solution
        p_ve = interpolate(p, self.Qe)
        p_exact_ve = interpolate(pe, self.Qe)
        e_p = abs(p_ve - p_exact_ve)

        # Define variational form of the error
        var_grad_e_p = inner(kappa * Grad(e_p, self.gdim), Grad(e_p, self.gdim))
        var_e_p = inner(e_p, e_p)

        # Assembling variational form of the error
        norm_grad_e_p = assemble(var_grad_e_p * dx(domain=self.mesh))
        norm_e_p = assemble(var_e_p * dx(domain=self.mesh))

        if test_params["error_format"] == "relative_error":
            norm_p = assemble( (float(tau / mu_f) * inner(kappa * Grad(p_exact_ve, self.gdim), Grad(p_exact_ve, self.gdim))
                               + float(beta) * inner(p_exact_ve, p_exact_ve)) * dx(domain=self.mesh))
        elif test_params["error_format"] == "absolute_error":
            norm_p = 1

        error_array[i] = (float(tau / mu_f) * norm_grad_e_p + float(beta) * norm_e_p)/ norm_p

        if  test_params["full_documentation"]:
            print '%--------------------'
            print '% Error in pressure:'
            print '%--------------------'
            print "\| grad e_p \|^2 = %.4e" % norm_grad_e_p
            print "\| e_p \|^2 = %.4e" % norm_e_p
            print "||| e_p |||^2 = %.4e" % error_array[i]

        return error_array, var_grad_e_p, var_e_p, norm_p


    def error_displacement(self, e_u_array, i, u, ue, mu, lmbda):
        """
        :param u: approximate solution
        :param ue: exact solution
        :param mu: first Lame coefficient
        :param lmbda: Lame coefficient
        :return: error-norm between u and ue
        """
        # Interpolate exact and approximate solution to the functional space of exact solution
        u_ve = interpolate(u, self.VVe)
        u_exact_ve = interpolate(ue, self.VVe)
        e = (u_ve - u_exact_ve)

        # Define variational form of the error
        var_strain_e_u = inner(strain_tensor(e, self.gdim), strain_tensor(e, self.gdim))
        var_div_e_u = inner(Div(e, self.gdim), Div(e, self.gdim))

        # Assembling variational form of the error
        strain_e_u = assemble(var_strain_e_u * dx(domain=self.mesh))
        div_e_u = assemble(var_div_e_u * dx(domain=self.mesh))

        if test_params["error_format"] == "relative_error":
            norm_u = assemble(2.0 * float(mu) * (inner(strain_tensor(u_exact_ve, self.gdim), strain_tensor(u_exact_ve, self.gdim)) \
                              + float(lmbda) * inner(Div(u_exact_ve, self.gdim), Div(u_exact_ve, self.gdim))) * dx(domain=self.mesh))
        elif test_params["error_format"] == "absolute_error":
            norm_u = 1

        if test_params["full_documentation"]:
            print '%------------------------'
            print '% Error in displacement:'
            print '%------------------------'
            print "\| eps(e_u) \|^2 = %.4e" % strain_e_u
            print "\| div(e_u) \|^2 = %.4e" % div_e_u
            print "||| e_u |||^2 = %.4e" % e_u_array[i]

        e_u_array[i] = (2.0 * float(mu) * strain_e_u + float(lmbda) * div_e_u) / norm_u

        return e_u_array, var_strain_e_u, var_div_e_u, norm_u

class BiotSolvers():

    # Init the BiotSolver with the mesh
    def __init__(self, mesh):
        self.mesh = mesh

    # Function with iterative coupling procedure
    def iterative_coupling(self, u, p, v, theta, bc_p, bc_u,
                                 error_control, e_p_array, e_u_array,
                                 phi_f, g, alpha, beta, mu, lmbda, L, mu_f, kappa, gdim, tau, n, ds,
                                 u_k, u_k1, p_k, p_k1,
                                 s_f_tilde, f, funcs,
                                 results_folder_name, test_params):
        ITER_NUM = test_params["iter_num"]

        for iter_counter in range(1, ITER_NUM):
            if test_params["full_documentation"]:
                print '%----------------------------'
                print '% Iteration # ', iter_counter
                print '%----------------------------'


            # ----------------------------------------------------------------------------------------------#
            # step 1: Constructing p (p_{k + 1}) using u_k
            # ----------------------------------------------------------------------------------------------#
            a_p_var = (L * inner(p, theta) + beta * inner(p, theta) \
                       + tau / mu_f * inner(kappa * Grad(p, gdim), Grad(theta, gdim))) * dx(domain=self.mesh)
            l_p_var = (inner(L * p_k - alpha * Div(u_k, gdim), theta) + inner(s_f_tilde, theta)
                       + phi_f * tau / mu_f * inner(kappa * g, Grad(theta, gdim))) * dx(domain=self.mesh)

            # Assemble the system to solve flow equation for pressure
            A_p, b_p = assemble_system(a_p_var, l_p_var, bc_p)
            solve(A_p, p_k1.vector(), b_p)

            # Calculate the error in the pressure term
            e_p_array, var_grad_e_p, var_e_p, norm_p = error_control.error_pressure(e_p_array, iter_counter, p_k1, funcs["p_e"], mu_f, kappa, tau, beta)

            #----------------------------------------------------------------------------------------------#
            # Constructin u (u_{k + 1}) using p_{k + 1}
            # ----------------------------------------------------------------------------------------------#
            a_u_var = (2.0 * mu * inner(strain_tensor(u, gdim), strain_tensor(v, gdim)) \
                       + lmbda * inner(Div(u, gdim), Div(v, gdim))) * dx(domain=self.mesh)
            l_u_var = (inner(f, v) - alpha * inner(Grad(p_k1, gdim), v)) * dx(domain=self.mesh) \
                      + inner(funcs["t_N"], v) * ds(4)

            # Assemble the system to solve mechanics equation for displacement
            A_u, b_u = assemble_system(a_u_var, l_u_var, bc_u)
            solve(A_u, u_k1.vector(), b_u)

            # Calculate the error in the displacement term
            e_u_array, var_strain_e_u, var_div_e_u, norm_u = error_control.error_displacement(e_u_array, iter_counter, u_k1, funcs["u_e"], mu, lmbda)

            # Update the iterations
            u_k.assign(u_k1)
            p_k.assign(p_k1)

        # Output the convergence of the pressure and displacement
        postprocess.output_errors_wrt_iterations(e_u_array, e_p_array, ITER_NUM)
        # Plot pressure
        postprocess.plot_function(p_k1, self.mesh, gdim,
                                  project_path + results_folder_name + 'pressure-n-%d-k-%d' % (n + 1, ITER_NUM))
        # Plot displacement
        postprocess.plot_vector_function(u_k1, self.mesh, gdim,
                                         project_path + results_folder_name + 'displacement-n-%d-k-%d' % (n + 1, ITER_NUM))

        return u_k1, p_k1, e_p_array, e_u_array

    # Function with full-implicit procedure
    def full_implicit_coupling(self, u, p, v, theta, bc_p, bc_u,
                                     alpha, beta, gdim, t_n1, tau,
                                     w_n1, s_f_tilde,
                                     mu_f, kappa, mu, lmbda,
                                     func_spaces, funcs, facet_function):
        # Collect BC for the mixed method
        bcs_p = [bc_p, bc_u]

        # Variational forms of the Biot RHS and LHS
        biot_lhs = (beta * inner(p, theta) + tau / mu_f * inner(kappa * Grad(p, gdim), Grad(theta, gdim))
                    + alpha * inner(theta, Div(u, gdim))
                    + 2.0 * mu * inner(strain_tensor(u, gdim), strain_tensor(v, gdim)) + lmbda * inner(Div(u, gdim), Div(v, gdim))
                    - alpha * inner(p, Div(v, gdim))) * dx(domain=self.mesh)
        biot_rhs = (inner(funcs["f"], v) + inner(s_f_tilde, theta)) * dx(domain=self.mesh)
        # Solve the Biot system
        solve(biot_lhs == biot_rhs, w_n1, bcs_p)
        u_n1, p_n1 = w_n1.split()

        return u_n1, p_n1

# Class of the Dirichlet BC
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

# Class for testing Biot problem
class TestBiot():
    # Init variable of TestBiot class
    def __init__(self, problem_data, domain_params, material_params, test_num):

        self.problem_data = problem_data    
        self.domain_params = domain_params
        self.material_params = material_params
        self.test_num = test_num

    def get_initial_pressure_and_displacement(self, mesh, ds, n,
                                              p, u, v, theta, bc_p, bc_u, p_k, u_k,
                                              f, t_N,
                                              phi_f, g, alpha, mu, lmbda):
        '''
        print '%--------------------'
        print '% Initial iteration:'
        print '%--------------------'
        '''
        # Define dim of problem
        gdim = self.domain_params['gdim']

        # Define the variational formulation for pressure
        a_p0 = inner(Grad(p, gdim), Grad(theta, gdim)) * dx(domain=mesh)
        l_p0 = phi_f * inner(g, Grad(theta, gdim)) * dx(domain=mesh)

        # Assemble the system for the initial pressure
        A_p0, b_p0 = assemble_system(a_p0, l_p0, bc_p)
        solve(A_p0, p_k.vector(), b_p0)

        """
        print "A_p0 = ", A_p0.array()
        print "b_p0 = ", b_p0.array()
        print "p0 = ",   p_k.vector().array()
        """
        # Define the variationalphi_0 formulation for displacement
        a_u0 = (2.0 * mu * inner(strain_tensor(u, gdim), strain_tensor(v, gdim)) \
                + lmbda * inner(Div(u, gdim), Div(v, gdim))) * dx(domain=mesh)
        l_u0 = (inner(f - alpha * Grad(p_k, gdim), v)) * dx(domain=mesh)

        # Assemble the system for initial displacement
        A_u0, b_u0 = assemble_system(a_u0, l_u0, bc_u)
        solve(A_u0, u_k.vector(), b_u0)

        """
        print "A_u0 = ", A_u0.array()
        print "b_u0 = ", b_u0.array()
        print "u0 = ", u_k.vector().array()
        """

        return p_k, u_k

    # Construction of the required functional space
    def functional_spaces(self, mesh, test_params):

        # Order of the space of the exact solutions
        v_exact_degree = test_params['u_approx_order'] + 2
        p_exact_degree = test_params['p_approx_order'] + 2
        
        # Space for the displacement
        V = VectorFunctionSpace(mesh, "CG", test_params['u_approx_order'])
        V_EL = FiniteElement("CG", mesh.ufl_cell(), test_params['u_approx_order'])

        # Space for the pressure
        Q = FunctionSpace(mesh, "CG", test_params['p_approx_order'])
        Q_EL = FiniteElement("CG", mesh.ufl_cell(), test_params['u_approx_order'])

        W = FunctionSpace(mesh, V_EL * Q_EL)

        # Exact space for the displacement
        V_exact = VectorFunctionSpace(mesh, "CG", v_exact_degree)
        Vex_EL = FiniteElement("CG", mesh.ufl_cell(), v_exact_degree)

        # Exact space for the pressure
        Q_exact = FunctionSpace(mesh, "CG", p_exact_degree)
        Qex_EL = FiniteElement("CG", mesh.ufl_cell(), p_exact_degree)
        W = FunctionSpace(mesh, V_EL * Q_EL)
        W_exact = Vex_EL * Qex_EL

        if self.domain_params['gdim'] == 1:
            H_div = FunctionSpace(mesh, 'CG', test_params['flux_approx_order'])
        else:
            H_div = FunctionSpace(mesh, 'RT', test_params['flux_approx_order'])

        func_spaces = dict(V = V, 
                           Q = Q,
                           W = W,
                           V_exact = V_exact,
                           Q_exact = Q_exact,
                           H_div = H_div,
                           W_exact = W_exact)
        return func_spaces

    def fetch_required_paramenters(self):

        E = float(self.material_params["E"])
        nu = float(self.material_params["nu"])
        alpha = float(self.material_params["alpha"])
        M = float(self.material_params["M"])
        c_f = float(self.material_params["c_f"])
        phi_0 = float(self.material_params["phi_0"])
        mu_f = float(self.material_params["mu_f"])
        lmbda = float(self.material_params["lmbda"])
        mu = float(self.material_params["mu"])

        # Convert matrix to UFL form
        if domain_params["gdim"] == 1:
            kappa = Expression(self.material_params["kappa"][0])
        else:
            kappa = as_matrix(self.material_params["kappa"])

        return E, nu, alpha, M, c_f, phi_0, mu_f, kappa, lmbda, mu

    def update_term_dependent_on_time(self, funcs, t_n):


        if self.test_num == 3:
            funcs['f'].t = t_n
            funcs['p_D_gamma1'].t = t_n
            funcs['s_f'].t = t_n

        else:
            funcs['u_e'].t = t_n
            funcs['p_e'].t = t_n
            funcs['f'].t = t_n
            funcs['s_f'].t = t_n

        return funcs

    def convert_extressions_to_functions(self, test_params):
        '''
        :return:
        '''
        if self.test_num == 1:

            # Exact pressure
            p_e = Expression(self.problem_data['p_expr'],
                             F = self.material_params['F'],
                             B = self.material_params['B'],
                             nu_u = self.material_params['nu_u'],
                             a = self.domain_params['l_x'],
                             alpha_1=self.material_params['alpha_n'][0],
                             alpha_2=self.material_params['alpha_n'][1],
                             alpha_3=self.material_params['alpha_n'][2],
                             alpha_4=self.material_params['alpha_n'][3],
                             alpha_5=self.material_params['alpha_n'][4],
                             c_f = self.material_params['c_f'],
                             t = 0.0, degree=3)
            # Initial condition for pressure
            p_0 = p_e
            # Dirichlet BC for pressure
            #p_D = Expression(self.problem_data['pD_expr'])
            # Flow source/sink
            s_f = Expression(self.problem_data['sf_expr'], degree=3)

            f = Expression((self.problem_data['f_expr'][0], self.problem_data['f_expr'][1]), degree=3)
            u_e = Expression((self.problem_data['u_expr'][0], self.problem_data['u_expr'][1]),
                              F = self.material_params['F'],
                              nu_u = self.material_params['nu_u'],
                              nu = self.material_params['nu'],
                              a = self.domain_params['l_x'],
                              mu = self.material_params['mu'],
                              alpha_1 = self.material_params['alpha_n'][0],
                              alpha_2 = self.material_params['alpha_n'][1],
                              alpha_3 = self.material_params['alpha_n'][2],
                              alpha_4 =self.material_params['alpha_n'][3],
                              alpha_5 =self.material_params['alpha_n'][4],
                              c_f = self.material_params['c_f'],
                              t = 0.0, degree=3)
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
                             t=0.0, degree=3)
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
                               t=0.0, degree=3)
            t_N = Expression((self.problem_data['t_N_expr'][0], self.problem_data['t_N_expr'][1]),
                             F=self.material_params['F'], degree=3)
            u_0 = u_e

            funcs = dict(u_e = u_e, u_e_0=u_e_0, u_e_1=u_e_1,
                    p_e = p_e,
                    u_0 = u_0,
                    p_0 = p_0,
                    f = f,
                    t_N = t_N,
                    s_f = s_f)

        elif self.test_num == 2:
            # Exact pressure and displacement
            p_e = Expression(self.problem_data['p_expr'], t=0.0, degree=3)
            u_e = Expression((self.problem_data['u_expr'][0], self.problem_data['u_expr'][1]), t=0, degree=3)

            # Initial condition for pressure and displacement
            p_0 = Expression(self.problem_data['p0_expr'], degree=3)
            u_0 = Expression((self.problem_data['u0_expr'][0], self.problem_data['u0_expr'][1]), degree=3)

            # Flow source/sink
            # sf_expr = '(1/M + c_f*phi_0) * x[0] * (1 - x[0]) * x[1] * (1 - x[1]) + 2 * t * (x[0] * (1 - x[0]) + x[1] * (1 - x[1])) + alpha * (1 + 2 * x[0])'
            s_f = Expression(self.problem_data['sf_expr'], c_f=self.material_params['c_f'],
                             phi_0=self.material_params['phi_0'], M=self.material_params['M'],
                             alpha=self.material_params['alpha'], t=0.0, degree=3)

            # Mechanical load
            # f_0_expr = '- mu * 4 * t - 2 * lmbda * t + alpha * x[1] * (1 - x[1]) * (1 - 2 * x[0]) * t'
            # f_1_expr = 'alpha * x[0] * (1 - x[0]) * (1 - 2 * x[1]) * t'
            f = Expression((self.problem_data['f_expr'][0], self.problem_data['f_expr'][1]),
                           mu = self.material_params['mu'],
                           lmbda = self.material_params['lmbda'],
                           alpha = self.material_params['alpha'],
                           t = 0.0, degree=3)
            t_N = Expression((self.problem_data['t_N_expr'][0], self.problem_data['t_N_expr'][1]), degree=3)

            funcs = dict(u_e=u_e,
                    p_e=p_e,
                    u_0=u_0,
                    p_0=p_0,
                    f=f,
                    s_f=s_f,
                    t_N=t_N)

        elif self.test_num == 3:
            # Dirichlet boundary conditions for pressure
            p_D_gamma1 = Expression(self.problem_data['p_D_expr'][0], t = 0.0)
            p_D_gamma = Expression(self.problem_data['p_D_expr'][1])

            # Initial condition for pressure
            p_0 = Expression(self.problem_data['p0_expr'])

            # Neumann boundary conditions for pressure
            t_N = Constant((self.problem_data['t_N_expr'][0], self.problem_data['t_N_expr'][1]))
            u_0 = Expression((self.problem_data['u0_expr'][0], self.problem_data['u0_expr'][1]))

            s_f = Expression(self.problem_data['sf_expr'])
            f = Expression(self.problem_data['f_expr'][0], self.problem_data['f_expr'][1])

            funcs = dict(p_D_gamma1 = p_D_gamma1,
                    p_D_gamma = p_D_gamma,
                    u_0 = u_0,
                    p_0 = p_0,
                    f = f,
                    s_f = s_f,
                    t_N=t_N)
        elif self.test_num == 4:
            # Exact pressure and displacement
            p_e = Expression(self.problem_data['p_expr'], t=0.0, degree=3)
            u_e = Expression((self.problem_data['u_expr'][0], self.problem_data['u_expr'][1]), t=0, degree=3)

            # Initial condition for pressure and for pressure
            p_0 = Expression(self.problem_data['p0_expr'], degree=3)
            u_0 = Expression((self.problem_data['u0_expr'][0], self.problem_data['u0_expr'][1]), degree=3)
            # Flow source/sink
            s_f = Expression(self.problem_data['sf_expr'], c_f=self.material_params['c_f'],
                             phi_0=self.material_params['phi_0'], M=self.material_params['M'],
                             alpha=self.material_params['alpha'], t=0.0, degree=3)
            # Mechanical load
            f = Expression((self.problem_data['f_expr'][0], self.problem_data['f_expr'][1]),
                           mu=self.material_params['mu'], lmbda=self.material_params['lmbda'],
                           alpha=self.material_params['alpha'], t=0.0, degree=3)
            t_N = Constant((self.problem_data['t_N_expr'][0], self.problem_data['t_N_expr'][1]))

            funcs = dict(p_e=p_e,
                         u_e=u_e,
                         u_0=u_0,
                         p_0=p_0,
                         f=f,
                         s_f=s_f, t_N=t_N)
        return funcs

    def test_biot(self, test_params, project_path):

        # Check if exist and create folder for the results
        results_folder_name = postprocess.construct_result_folder_name(self.domain_params['gdim'], self.test_num, test_params)
        postprocess.create_results_folder(results_folder_name)

        # Generate mesh and boundary parts
        if self.domain_params['gdim'] == 2:
            facet_function = problem.geometry_2d(self.domain_params, self.test_num, test_params["mesh_resolution"])
        elif self.domain_params['gdim'] == 3:
            facet_function = problem.geometry_3d(self.domain_params, self.test_num)
        # Get the mesh from facet_function
        mesh = facet_function.mesh()
        if self.domain_params['gdim'] == 2:
            postprocess.plot_mesh(mesh, project_path + results_folder_name + 'initial-mesh')
        elif self.domain_params['gdim'] == 3:
            postprocess.plot_mesh_3d(mesh, project_path + results_folder_name + 'initial-mesh')

        # Calculate the estimate for Freidrichs constant based on the domain
        C_FD = problem.calculate_CF_of_domain(domain_params)

        # Construct required functions 
        funcs = self.convert_extressions_to_functions(test_params)

        # Construct reuired functional spaces
        func_spaces = self.functional_spaces(mesh, test_params)

        # Set Dirichlet BC for the problem,
        bc_p, bc_u = self.dirichlet_boundary_conditions(func_spaces, funcs, facet_function, material_params)

        postprocess.output_problem_characteristics(test_num, mesh, problem_data, domain_params, material_params, test_params, func_spaces)

        # Set start time
        t0_problem = time.clock()
        # Solve Biot system
        self.solve_biot(func_spaces, funcs, facet_function, mesh, C_FD, bc_p, bc_u,
                        test_params,
                        project_path, results_folder_name)
        # Set the final time
        tfinal_problem = time.clock()
        print "total time = %d mins" %(float(tfinal_problem - t0_problem) / 60.0)

        interactive()

    def solve_biot(self, func_spaces, funcs, facet_function, mesh, C_FD, bc_p, bc_u,
                   test_params, 
                   project_path, results_folder_name):

        # Fetch required parameters
        E, nu, alpha, M, c_f, phi_0, mu_f, kappa, lmbda, mu = self.fetch_required_paramenters()
        beta = 1 / M + c_f * phi_0 #
        if self.test_num == 1: # Mandel's problem
            g = Constant((0, 0.0))
        else:
            g = Constant((0, -9.8))
        phi_f = phi_0

        print 'alpha = ', float(alpha)
        print 'beta = ', float(beta)
        print "M = ", float(M)
        print "E = ", float(E)
        print "nu = ", float(nu)
        print "mu_f = ", float(mu_f)

        # Contractive algorith parameters
        L = Constant(alpha**2 / (2 * lmbda))
        q = Constant(L / (L + beta))

        print 'L = ', float(L)
        print 'q = ', float(q)

        # Define the value of the time-step
        tau = float(self.problem_data['t_T'] / test_params['time_steps'])

        # Initialize times
        t_n = 0     # t_n
        t_n1 = tau  # t_{n + 1}

        # Initialial time counter
        n = 0
        
        # Geometry parameters
        ds = Measure("ds", subdomain_data=facet_function)
        normal = FacetNormal(mesh)
        gdim = self.domain_params['gdim']

        # Define function p_{n-1} and u_{n-1}, which is p_0 and u_0 on the initial step
        u_n = interpolate(funcs['u_0'], func_spaces['V'])
        p_n = interpolate(funcs['p_0'], func_spaces['Q'])

        # Define the source function on the n-th time step
        s_f_tilde = tau * funcs["s_f"] + beta * p_n + alpha * Div(u_n, gdim)

        # Define the error control class
        error_control = ErrorControl(mesh, func_spaces['V_exact'], func_spaces['Q_exact'],
                                     self.domain_params["gdim"])
        e_u_array = error_control.get_error_array(test_params["iter_num"])
        e_p_array = error_control.get_error_array(test_params["iter_num"])

        # Define the boit solver class
        biot_solver = BiotSolvers(mesh)

        # Loop over time-steps
        while t_n1 < self.problem_data['t_T'] + DOLFIN_EPS:

            print '%----------------------------------------------------'
            print ' Time step (t_%i, t_%i) = (%.3f, %.3f):' %(n, n+1, t_n, t_n1)
            print '%----------------------------------------------------'

            # Update the fource term
            funcs = self.update_term_dependent_on_time(funcs, t_n1)

            # Update Dirichlet BC if the given data is dependent on time
            bc_p, bc_u = self.dirichlet_boundary_conditions(func_spaces, funcs, facet_function, material_params)

            postprocess.plot_function(interpolate(funcs["p_e"], func_spaces["Q_exact"]),
                                      mesh, gdim,
                                      project_path + results_folder_name + 'pressure-exact-n-%d' % (n + 1))
            postprocess.plot_vector_function(interpolate(funcs["u_e"], func_spaces["V_exact"]),
                                      mesh, gdim,
                                      project_path + results_folder_name + 'displacement-exact-n-%d' % (n + 1))

            if test_params['coupling_approach'] == 'fully-implicit':

                u, p = TrialFunctions(func_spaces["W"])
                v, theta = TestFunctions(func_spaces["W"])
                w_n1 = Function(func_spaces["W"])

                u_n1, p_n1 = biot_solver.full_implicit_coupling(u, p, v, theta, bc_p, bc_u,
                                                                alpha, beta, gdim, t_n1, tau,
                                                                w_n1, s_f_tilde,
                                                                mu_f, kappa, mu, lmbda,
                                                                func_spaces, funcs, facet_function)
                # Calculate the error in the pressure term
                e_p_array = error_control.error_pressure(e_p_array, 0, p_n1, funcs["p_e"], mu_f, kappa, tau, beta)
                e_p_array = error_control.error_displacement(e_p_array, 0, u_n1, funcs["u_e"], mu, lmbda)

            else:
                # Define trial fucntions for mechanics and flow equations
                u = TrialFunction(func_spaces['V'])
                p = TrialFunction(func_spaces['Q'])

                # Define test functions for mechanics and flow equations
                v = TestFunction(func_spaces['V'])
                theta = TestFunction(func_spaces['Q'])

                # Define function p_{n, k},p_{n, k+1}  and u_{n, k}, u_{n, k+1}
                u_k = Function(func_spaces['V'])
                u_k1 = Function(func_spaces['V'])
                p_k = Function(func_spaces['Q'])
                p_k1 = Function(func_spaces['Q'])

                # Set initial values of pressure and displacement
                p_k = interpolate(Expression("0.0"), func_spaces['Q'])
                u_k = interpolate(Expression(("0.0", "0.0")), func_spaces['V'])

                '''
                p_k, u_k = self.get_initial_pressure_and_displacement(mesh, ds, normal,
                                                                      p, u, v, theta, bc_p, bc_u,
                                                                      p_k, u_k,
                                                                      funcs['f'], funcs['t_N'],
                                                                      phi_f, g, alpha, mu, lmbda)
                '''
                """
                u_k = u_n
                p_k = p_n
                """

                # Plot pressure
                postprocess.plot_function(p_k, mesh, gdim,
                                          project_path + results_folder_name + 'pressure-n-%d-k-%d' % (n+1, 0))
                postprocess.plot_vector_function(u_k, mesh, gdim,
                                                 project_path + results_folder_name + 'displacement-n-%d-k-%d' % (n+1, 0))

                # Calculate the error in the pressure term
                e_p_array, var_grad_e_p, var_e_p, norm_p = error_control.error_pressure(e_p_array, 0, p_k, funcs["p_e"], mu_f, kappa, tau, beta)
                e_u_array, var_strain_e_u, var_div_e_u, norm_u = error_control.error_displacement(e_u_array, 0, u_k, funcs["u_e"], mu, lmbda)

                u_n1, p_n1, e_p_array, e_u_array = biot_solver.iterative_coupling(u, p, v, theta, bc_p, bc_u,
                                                                                error_control, e_p_array, e_u_array,
                                                                                phi_f, g, alpha, beta, mu, lmbda, L, mu_f, kappa, gdim, tau, n, ds,
                                                                                u_k, u_k1, p_k, p_k1,
                                                                                s_f_tilde, funcs["f"], funcs,
                                                                                results_folder_name, test_params)

            # Plot pressure
            postprocess.plot_function(p_n1, mesh, gdim,
                                      project_path + results_folder_name + 'pressure-n-%d' % (n + 1))
            # Plot displacement
            postprocess.plot_vector_function(u_n1, mesh, gdim,
                                             project_path + results_folder_name + 'displacement-n-%d' % (n + 1))

            # Update time layer
            t_n1 += tau
            t_n += tau
            n += 1

            # Update functions
            u_n.assign(u_n1)  # u_{n} <- u_{n + 1}
            p_n.assign(p_n1)  # p_{n} <- p_{n + 1}

            # Set Dirichlet BC for pressure and displpacement

    def dirichlet_boundary_conditions(self, func_spaces, funcs, facet_function, material_params):

        # left: 1
        # right: 2
        # bottom: 3
        # top: 4

        if self.test_num == 1:
            # funcs['u_D'], funcs['p_D'],
            # Set Dirichlet BC for pressure on the right(2) boundary of the domain
            bcs_p = [DirichletBC(func_spaces['Q'], Constant(0.0), facet_function, 2)]

            # Set Dirichlet BC for displacement on the left(1) and bottom(3) boundary of the domain
            bcs_u = [DirichletBC(func_spaces['V'].sub(0), Constant(0.0) , facet_function, 1), # u_x = 0
                     DirichletBC(func_spaces['V'].sub(1), Constant(0.0), facet_function, 3), # u_y = 0
                     DirichletBC(func_spaces['V'].sub(1), funcs["u_e_1"], facet_function, 4)]  # u_y = U_2
            """
            bcs_u = [DirichletBC(func_spaces['V'].sub(0), Constant(0.0), facet_function, 1),
                     DirichletBC(func_spaces['V'].sub(1), Constant(0.0), facet_function, 3),
                     DirichletBC(func_spaces['V'].sub(1), Constant(2*material_params["F"]), facet_function, 4)]
            """
        elif self.test_num == 2:
            # Set Dirich BC for pressure and displacement on the whole boundary
            bcs_p = DirichletBC(func_spaces['Q'], funcs['p_e'], DirichletBoundary())
            bcs_u = DirichletBC(func_spaces['V'], funcs['u_e'], DirichletBoundary())

        elif self.test_num == 3:
            bcs_p = [DirichletBC(func_spaces['Q'], funcs['p_e'], facet_function, 1),
                     DirichletBC(func_spaces['Q'], funcs['p_e'], facet_function, 0)]
            bcs_u = []
        elif self.test_num == 4:
            # Set Dirich BC for pressure and displacement on the whole boundary
            bcs_p = DirichletBC(func_spaces['Q'], funcs['p_e'], DirichletBoundary())
            bcs_u = DirichletBC(func_spaces['V'], funcs['u_e'], DirichletBoundary())

        return bcs_p, bcs_u


if __name__ == '__main__':

    # Defnie the file path and project path
    project_path = postprocess.get_project_path()
        
    # Types of domain
    rectangular_domain_tag = "rectangular-domain"
    rectangular_domain_with_obsticle_tag = "rectangular-with-obstacle"

    # Decoupling
    fully_implicit_coupling = 'fully-implicit'
    explicit_coupling = 'explicit'
    iterative_implicit_coupling = 'iterative'

    # Error format
    relative_error = "relative_error"
    absolute_error = "absolute_error"

    # Pressure recovery method
    CG_method = 'CG-method'
    mixed_method = 'mixed-method'

    # Init problem parameters
    test_params = dict(u_approx_order = 1,        # Functional spaces parameters for displacement
                       p_approx_order = 1,        # Functional spaces parameters for pressure
                       flux_approx_order = 2,
                       iter_accuracy = 1e-4,    # Required accuracy at each interation cycle
                       time_steps = 40,         # Number of time steps on the interval [0, t_T]
                       mesh_resolution = 20,    # Lever of refiniment of initial mesh
                       iter_num = 5,
                       coupling_approach = iterative_implicit_coupling,
                       pressure_recovery_method = CG_method,
                       full_documentation=False,
                       error_format=absolute_error)
    # list of tests
    tests = {1: tests.mandels_problem_2d_t,
             2: tests.simple_example_2d_t,
             3: tests.kolesov_example_2d_t,
             4: tests.manuel_example_2d_t}

    # Set the number of the test and call for the problem data
    test_num = 1
    problem_data, domain_params, material_params = tests[test_num]()

    test = TestBiot(problem_data, domain_params, material_params, test_num)
    test.test_biot(test_params, project_path)

