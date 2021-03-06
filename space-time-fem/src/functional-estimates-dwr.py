__author__ = "Svetlana Matculevich <smatculevich@ricam.oeaw.ac.at>"

from dolfin import *
from dolfin.cpp.function import near
from mshr import generate_mesh, Box, Cone, Rectangle
from dolfin.cpp.mesh import UnitSquareMesh, Mesh, UnitIntervalMesh, UnitCubeMesh, Point, RectangleMesh, FacetFunction, \
    CellFunction, SubDomain
import time
import tests, problem, estimates, postprocess
from dolfin.cpp.common import set_log_level
from problem import Grad, NablaGrad, dirichlet_bc_marker, neumann_bc_marker, init_bc_marker, do_nothing_bc_marker


CRITICAL  = 50 # errors that may lead to data corruption and suchlike
ERROR     = 40 # things that go boom
WARNING   = 30 # things that may go boom later
INFO      = 20 # information of general interest
PROGRESS  = 16 # what's happening (broadly)
TRACE     = 13 # what's happening (in detail)
DBG       = 10  # sundry
set_log_level(20)

parameters["plotting_backend"] = "matplotlib"
parameters ["form_compiler"]["cpp_optimize"] = True
parameters['allow_extrapolation'] = True
#parameters["refinement_algorithm"] = "plaza_with_different_faces"
#parameters["refinement_algorithm"] = "plaza_with_parent_faces"

expression_degree = 4

# Define Dirichlet boundary condition
def boundary(x, on_boundary):
    return on_boundary

# Class of the Dirichlet BC
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


#----------------------------------------------------------------------------------------------------------------------#
# Class TestEstimates
#----------------------------------------------------------------------------------------------------------------------#
class TestEstimates():

    def __init__(self, test_num, u_expr, grad_u_expr, f_expr, 
                       A_expr, min_eig_A, lmbd_expr, a_expr, c_H, eps,
                       u0_expr, uD_expr, uN_expr, T, domain, dim):

        # Initialize the problem data
        self.u_expr = u_expr
        self.grad_u_expr = grad_u_expr
        self.f_expr = f_expr
        self.A_expr = A_expr
        self.min_eig_A = min_eig_A
        self.lmbd_expr = lmbd_expr
        self.a_expr = a_expr
        self.c_H = c_H
        self.eps = eps
        self.T = T
        self.u0_expr = u0_expr
        self.uD_expr = uD_expr
        self.uN_expr = uN_expr

        self.boundary = boundary
        self.domain = domain
        self.dim = dim

        self.test_num = test_num

    # Convert given data expressions form string format to the UFL Expressions
    def convert_problem_data_to_expressions(self, mesh, problem_params):
        """
        :param mesh: mesh-discretization of the domain
        :param dim: dimension of the problem
        :return: given data as functions
                 f, A, invA, adjA, lmbd, a, uD, u_e, grad_u_e
        """

        # Define diffusion matrix A if it is constant over the whole domain
        if problem_params['material_tag'] == "material-constant":
            # Convert matrix to UFLMatrix
            if dim == 2:
                A = Constant(self.A_expr)
            else:
                A = as_matrix(self.A_expr)
        else:
            A = Expression(self.A_expr, degree=expression_degree)

        if dim == 2:
            if problem_params['material_tag'] == "material-constant":
                invA = 1/A
            else:
                invA = Expression("1 / (" + self.A_expr + ")", degree=expression_degree)

        elif dim == 3: 
            #invA = [[1/A[0][0], 1/A[0][1], 0], [1/A[1][0], 1/A[1][1], 0], [0, 0, 0]]
            invA = [[1 / A[0][0], 0, 0], [0, 1 / A[1][1], 0], [0, 0, 0]]
        # Change that later
        adjA = invA

        # Construct expression for the convection/advection vector field
        if dim == 2:
            #a = Constant(self.a_expr[0])
            #a = Expression((self.a_expr, "0.0"), degree=expression_degree)
            a = Expression((self.a_expr), degree=expression_degree)
        elif dim == 3:
            a = Expression((self.a_expr[0], self.a_expr[1], "0.0"), degree=expression_degree)

        # Construct expression for the exact solution expression
        if problem_params['solution_tag'] == "reference-solution":
            u_e = "0"
            grad_u_e = []
        elif problem_params['solution_tag'] == "predefined-solution":
            u_e = Expression(self.u_expr, degree=expression_degree)
            if self.dim == 2:
                grad_u_e = Expression((self.grad_u_expr[0], "0.0"), degree=expression_degree)
            elif self.dim == 3:
                grad_u_e = Expression((self.grad_u_expr[0], self.grad_u_expr[1], "0.0"), degree=expression_degree)

        # Define initial state
        u0 = Expression(self.u0_expr, degree=expression_degree)
        uN = Expression(self.uN_expr, degree=expression_degree) # degree=problem_params["v_approx_order"]
        uD = u_e
        lmbd = Expression(self.lmbd_expr, degree=expression_degree)

        # Define right-hand side
        f = Expression(self.f_expr, degree=expression_degree)

        return f, A, invA, adjA, lmbd, a, uD, uN, u0, u_e, grad_u_e


    # Function to generate the mesh
    def construct_mesh(self, nx0, nx1, nx2, nt, res, project_path):
        """
        :param nx0, nx1, nx2: mesh-sizes wrt to X, Y, Z directions
        :param res: mesh resolution if defined with mesh_generator function
        :param project_path: path of the project to load the mesh
        :return: mesh:discretization of the domain
        """
        # --------------------------------------------------------------------------------------------------------------#
        if self.domain == "unit-domain":

            if self.dim == 2:
                # Create mesh on the unit square
                mesh = UnitSquareMesh(nx0, nt)

            elif self.dim == 3:
                # Create mesh on the unit cube
                mesh = UnitCubeMesh(nx0, nx1, nt)

        elif self.domain == "10byT-domain":

            if self.dim == 2:
                # Create mesh on the unit square
                point_1 = Point(0.0, 0.0)
                point_2 = Point(10.0, self.T)
                mesh = RectangleMesh(point_1, point_2, nx0, nt, 'left')

        elif self.domain == "long-rectangle-domain":
            
            if self.dim == 2:
                
                point_1 = Point(0.0, 0.0)
                point_2 = Point(1.0, 3.0)
                mesh = RectangleMesh(point_1, point_2, nx0, nt, 'left')

        elif self.domain == "rectangular-domain-1x2":
          
            if self.dim == 2:
                
                point_1 = Point(0.0, 0.0)
                point_2 = Point(1.0, 2.0)
                mesh = RectangleMesh(point_1, point_2, nx0, nt, 'left')

        return mesh
    
    def construct_boundary_faces(self, mesh):

        # Get boundary facets function
        facet_funcs = FacetFunction('size_t', mesh)
        facet_funcs.set_all(0)

        # Dirichlet part of the boudary is marked by 0
        # dirichlet_bc_marker = 1

        # By default left and right is Dirichlet BC
        left = AutoSubDomain(lambda x: near(x[0], 0.0))  # x_min
        right = AutoSubDomain(lambda x: near(x[0], 1.0))  # x_max

        left.mark(facet_funcs, dirichlet_bc_marker)
        right.mark(facet_funcs, dirichlet_bc_marker)

        if dim == 3:
            # By default left and right is Dirichlet BC
            top = AutoSubDomain(lambda x: near(x[1], 0.0))  # x_min
            bottom = AutoSubDomain(lambda x: near(x[1], 1.0))  # x_max

            top.mark(facet_funcs, dirichlet_bc_marker)
            bottom.mark(facet_funcs, dirichlet_bc_marker)

        #neumann_bc_marker = dirichlet_bc_marker + 1

        if self.test_num == 24:
            left = AutoSubDomain(lambda x: near(x[0], 0.0))     # x_min
            right = AutoSubDomain(lambda x: near(x[0], 10.0))    # x_max

            left.mark(facet_funcs, neumann_bc_marker)
            right.mark(facet_funcs, neumann_bc_marker)

        if self.test_num == 37:
            right = AutoSubDomain(lambda x: near(x[0], 1.0))    # x_max
            right.mark(facet_funcs, neumann_bc_marker)

        if self.test_num == 38 or test_num == 39 or test_num == 40 or test_num == 41 or test_num == 42:
            left = AutoSubDomain(lambda x: near(x[0], 0.0))  # x_min
            right = AutoSubDomain(lambda x: near(x[0], 1.0))    # x_max

            left.mark(facet_funcs, neumann_bc_marker)
            right.mark(facet_funcs, neumann_bc_marker)

        # The boundary with IC
        #init_bc_marker = neumann_bc_marker + 1 # = 3
        sigma_0 = AutoSubDomain(lambda x: near(x[1], 0.0))
        sigma_0.mark(facet_funcs, init_bc_marker)

        # The top of the space-time cylinder
        #do_nothing_bc_marker = init_bc_marker + 1 # = 4
        sigma_T = AutoSubDomain(lambda x: near(x[1], self.T))
        sigma_T.mark(facet_funcs, do_nothing_bc_marker)

        ds = Measure('ds')[facet_funcs]
        dS = Measure('dS')[facet_funcs]

        return facet_funcs, ds, dS

    def execute_refinement_strategy(self, problem_params, mesh, e_distr, md_distr,
                                            project_path, results_folder, init_data_info):
        # error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_max_k, h_min_k, t_T, k, nt,

        # Define the refinement strategy
        if problem_params['refinement_tag'] == "adaptive":
            # Define the distribution (error or majorant) upon which the refinement is basedy^
            if problem_params['refinement_criteria_tag'] == 'error':
                distr = e_distr
            elif problem_params['refinement_criteria_tag'] == 'majorant':
                distr = md_distr

            if problem_params['marking_tag'] == 'average':
                marking = estimates.averaged_marking(mesh, distr)
                refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params[
                    'marking_tag'] + '-marking'

            elif problem_params['marking_tag'] == 'predefined':
                marking = estimates.predefined_amount_of_elements_marking(mesh, distr,
                                                                          problem_params['percentage_value'])
                refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params[
                    'marking_tag'] + '-marking' + '-proc-%d' % (problem_params['percentage_value'] * 100)

            elif problem_params['marking_tag'] == 'bulk':
                marking = estimates.bulk_marking(mesh, distr, problem_params['percentage_value'])
                refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params[
                    'marking_tag'] + '-marking' + '-bulk-%d' % (problem_params['percentage_value'] * 100)
            print("% ----------------")
            print("% Refine the mesh ")
            print("% ----------------")

            # Refine mesh based on the marking criteria
            mesh = refine(mesh, marking)
            print('mesh parameters: num_cells = %d, num_vertices = %d' % (mesh.num_cells(), mesh.num_vertices()))

            # Define the name of the mat-file
            mat_file_tag = project_path + results_folder + 'results-on-mesh-based-' + refinement_tag + init_data_info

        elif problem_params['refinement_tag'] == "uniform":

            cell_markers = CellFunction("bool", mesh)
            cell_markers.set_all(True)
            mesh = refine(mesh, cell_markers)
            refinement_tag = problem_params['refinement_tag']

            # Define the name of the matfile
            mat_file_tag = project_path + results_folder + 'results-on-mesh-uniform' + init_data_info

        return mesh, mat_file_tag

    def test_estimates(self, test_params):

        # Define the project path
        # --------------------------------------------------------------------------------------------------------------#
        project_path = postprocess.get_project_path()

        # Check if it does exist and create folder for the results
        # --------------------------------------------------------------------------------------------------------------#
        results_folder_name = postprocess.construct_result_folder_name(self.dim, test_num, test_params)
        postprocess.create_results_folder(results_folder_name, test_params)

        # --------------------------------------------------------------------------------------------------------------#
        # Define problem data
        # --------------------------------------------------------------------------------------------------------------#

        # Define the mesh based on the domain geometry
        # --------------------------------------------------------------------------------------------------------------#
        mesh = self.construct_mesh(test_params["nx0"], test_params["nx1"], test_params["nx2"], test_params["nt"], test_params["res"], project_path)

        if test_params["PLOT"] == True:
            # Plot and save the initial mesh
            file_name = project_path + results_folder_name + 'initial-mesh'
            if self.dim == 2:   postprocess.plot_mesh(mesh, file_name)
            elif self.dim == 3: postprocess.plot_mesh_3d(mesh, file_name)

        boundary_facets, ds, dS = self.construct_boundary_faces(mesh)

        # Calculate the estimate for Freidrichs constant based on the domain
        C_FD = problem.calculate_CF_of_domain(domain, self.dim-1)
        C_Ntr = problem.calculate_trace_constant(self.test_num)

        # Define functional spaces based on the mesh
        #--------------------------------------------------------------------------------------------------------------#
        V, VV, V_exact, VV_exact, H_div = \
            problem.functional_spaces(mesh, test_params, self.dim)

        
        # Define problem data functions
        #--------------------------------------------------------------------------------------------------------------#
        f, A, invA, adjA, lmbd, a, uD, uN, u0, u_e, grad_u_e = self.convert_problem_data_to_expressions(mesh, test_params)

        # Majorant parameters
        delta = 1
        # gamma = 1.5
        # gamma = 2
        gamma = 1

        # Output to console the problem data
        problem.output_problem_characteristics(self.test_num, self.u_expr, self.grad_u_expr, self.f_expr,
                                               self.domain, self.dim, mesh.num_cells(), mesh.num_vertices(),
                                               test_params["v_approx_order"])

        t0_problem = time.clock()
        self.solve_problem_nd_t(mesh, boundary_facets, ds, dS,
                           V, VV, V_exact, VV_exact, H_div, 
                           u_e, grad_u_e, f, A, invA, adjA, self.min_eig_A, lmbd, a,
                           u0, uD, uN, C_FD, C_Ntr,
                           delta, gamma,
                           project_path, results_folder_name,
                           test_params)
        t1_problem = time.clock()

        print("time = %d s" % (t1_problem - t0_problem))

    #-------------------------------------------------------------------------------#
    def solve_problem_nd_t(self, mesh, boundary_facets, ds, dS,
                           V, VV, V_exact, VV_exact, H_div,
                           u_e, grad_ue, f, A, invA, adjA, min_eig_A, lmbd, a,
                           u0, uD, uN, C_FD, C_Ntr,
                           delta, gamma,
                           project_path, results_folder, test_params):

        # Define the number of refinements
        ref_num = test_params['number_of_refinements']


        ref_num_array = postprocess.allocate_array(ref_num)
        h_max_array = postprocess.allocate_array(ref_num)
        h_min_array = postprocess.allocate_array(ref_num)
        DOFs        = postprocess.allocate_array(ref_num)
        numcells_array = postprocess.allocate_array(ref_num)
        numverts_array = postprocess.allocate_array(ref_num)

        e_H1_array = postprocess.allocate_array(ref_num)
        e_L2_array = postprocess.allocate_array(ref_num)
        e_total_array = postprocess.allocate_array(ref_num)

        maj_array = postprocess.allocate_array(ref_num)
        maj_II_array = postprocess.allocate_array(ref_num)
        min_array = postprocess.allocate_array(ref_num)

        i_eff_maj_array = postprocess.allocate_array(ref_num)
        i_eff_maj_II_array = postprocess.allocate_array(ref_num)
        i_eff_min_array = postprocess.allocate_array(ref_num)

        rel_error_k = postprocess.allocate_array(ref_num)
        rel_majorant_k = postprocess.allocate_array(ref_num)

        e_L_array = postprocess.allocate_array(ref_num)
        e_id_array = postprocess.allocate_array(ref_num)

        # Initialize the variables before the loop
        maj = 10e8
        i = 0

        while i <= ref_num - 1:

            print(" ")
            print("%-----------------------------------------------------------------------------------------------------------%")
            print " Refinement cycle # %d : DOFs = %d" %(i, len(V.dofmap().dofs()))
            print("%-----------------------------------------------------------------------------------------------------------%")
            print(" ")

            # Compute approximate solution, its error and majorant
            u, delta_stab = problem.solve_convection_reaction_diffusion(V, c_H, eps, f, A, min_eig_A, lmbd,
                                                            problem.interpolate_vector_function(a, dim, V_exact, VV_exact),
                                                            u0, uD, uN,
                                                            mesh, boundary_facets,
                                                            self.dim, self.test_num, test_params)
            #u = problem.solve_parabolic_problem(V, c_H, eps, f, A, u0, uD, mesh, boundary_facets, dim, test_num)
            # In case of explicitely unknown exact solution, calculate it with high order order polynomials
            if test_params['solution_tag'] == "reference-solution":
                # Compute reference solution with high order polynomials
                u_e, delta_stab_e = problem.solve_convection_reaction_diffusion(V_exact, c_H, eps, f, A, min_eig_A, lmbd,
                                                                  problem.interpolate_vector_function(a, dim, V_exact, VV_exact),
                                                                  u0, uD, uN,
                                                                  mesh, boundary_facets,
                                                                  self.dim, self.test_num, test_params)
                #u_e = problem.solve_parabolic_problem(V_exact, c_H, eps, f, u0, uD, mesh, boundary_facets, dim, test_num)
            #'''
            if test_params["PLOT"] == True and ref_num <= 3:
                # Plot approximate solution
                postprocess.plot_function_3d(mesh, u, project_path + results_folder + 'u-%d' % i)
                #postprocess.plot_function_3d(mesh, u_e, project_path + results_folder + 'ue-%d' % i)

                #postprocess.plot_function(u_e, mesh, dim, project_path + results_folder + 'u-ref-%d' % i)
                #postprocess.plot_function_3d(mesh, f, project_path + results_folder + 'ue-ref-%d' % i)
                #postprocess.plot_function(u, mesh, dim, project_path + results_folder + 'u-%d' % i)
            #'''
            # Calculate error
            e, val_e, val_grad_e, val_delta_e, e_T, val_e_L, val_e_id, \
            var_e, var_delta_e, var_grad_e = estimates.error_norm(u, u_e, A, lmbd,
                                                                  problem.interpolate_vector_function(a, dim, V_exact, VV_exact),
                                                                  f, c_H,
                                                                  test_params["v_approx_order"],
                                                                  mesh, T, dim, V, V_exact)

            # Define the error for the majorant
            error = (2 - delta) * val_grad_e + 2 * val_delta_e + c_H * e_T
            error_II = (2 - delta) * val_grad_e + 2 * val_delta_e + (1 - 1/gamma) * c_H * e_T

            if test_params["error_estimates"] == True:

                # Calculate majorant
                #y_ = project(A * NablaGrad(u_e, dim), V)
                #y = project(A * NablaGrad(u, dim), H_div)

                #y = project(A * NablaGrad(interpolate(u_e, V_exact), dim), H_div)
                y = project(A * NablaGrad(u, dim), H_div)

                """
                if test_params['solution_tag'] == 'predefined-solution':
                    y = project(A * NablaGrad(u, dim), H_div)
                elif test_params['solution_tag'] == 'reference-solution':
                    y = project(A * NablaGrad(interpolate(u_e, V_exact), dim), H_div)
                """

                #print ' y - y_ = ', assemble(inner(y - y_, y - y_) * dx)
                #y = interpolate(problem.Grad(u, dim), H_div)
                MAJORANT_OPTIMIZE = 1

                #y = project(problem.Grad(u_e, dim), H_div)
                #MAJORANT_OPTIMIZE = 1

                # Calculate majorant
                maj, y, beta, md, mf, var_m_d, var_m_f_w_opt, \
                majorant_reconstruction_time, majorant_minimization_time  = \
                    estimates.majorant_nd(u, V_exact, y, H_div,
                                          f, A, invA, min_eig_A, c_H, eps, lmbd, a,
                                          u0, error, mesh, C_FD, dim, test_params)

                i_eff_maj = sqrt(maj / error)

                maj_array[i] = maj
                i_eff_maj_array[i] = i_eff_maj

                #w = project(u_e - interpolate(u, W), W)
                # Calculate majorant
                #maj_II, beta, majorant_reconstruction_time, majorant_minimization_time = \
                #    estimates.majorant_II_nd(u, V_exact, w, W, y, f, u0, error_II, gamma, T, mesh, C_FD, dim, v_deg, MAJORANT_OPTIMIZE)
                #i_eff_maj_II = sqrt(maj_II / error)

                # Construct error and majorant distribution
                #e_distr, m_distr, maj_distr = estimates.majorant_distribution_using_DG0(e, rd, rf, A, invA, beta, C_FD, mesh, dim)

                # , mf_distr, maj_distr, ed, md = \
                ed_distr, md_distr, maj_distr = \
                     estimates.error_majorant_distribution_nd(mesh, dim, V,
                                                              var_grad_e, var_delta_e,
                                                              var_m_d, var_m_f_w_opt, beta, C_FD)

                #min, phi = estimates.minorant(u, mesh, V_exact, u_0, boundary, f, dim, error)
                #i_eff_min = sqrt(min / error)
                #estimates.output_result_minorant(min, min/error)
                min = 0.0
                i_eff_min = 0.0

                #residual_CG0_distr, e_DWR_CG0_distr, J_e_CG0_distr = estimates.get_indicators_CG0(e, sqrt(error), V, V_exact, f, u0, boundary, u, u_e,
                #                                                                                  mesh, dim)
                #eta_distr_, E_DWR_distr_, J_e_distr_, m_d_distr_, m_df_distr_ = \
                #    majorants.compare_error_indicators(mesh, V, V_exact, V_exact, f, u0, boundary, u, interpolate(u_e, V_exact), y, beta, test_num)

            # Update the arrays of data with respect to the refinement cycle
            ref_num_array[i] = i

            e_H1_array[i] = sqrt(val_grad_e)
            e_L2_array[i] = sqrt(val_e)
            e_total_array[i] = sqrt(error)
            e_L_array[i] = sqrt(val_e_L)
            e_id_array[i] = sqrt(val_e_id)

            maj_array[i] = sqrt(maj)
            i_eff_maj_array[i] = sqrt(maj / error)

            maj_array[i] = sqrt(maj)
            i_eff_min_array[i] = sqrt(min / error)

            h_max_array[i] = mesh.hmax()
            h_min_array[i] = mesh.hmin()

            DOFs[i] = len(V.dofmap().dofs())
            num_cells = mesh.num_cells()
            num_verts = mesh.num_vertices()

            numcells_array[i] = num_cells
            numverts_array[i] = num_verts

            # Output the results of upper bound reconstructions
            estimates.output_result_error_and_majorant(sqrt(error), sqrt(maj), sqrt(maj / error))

            # Construct the tag with problem information
            results_info = postprocess.construct_result_tag(test_num, i, test_params["nx0"], test_params["nx1"],
                                                            test_params["nx2"],
                                                            num_cells, num_verts)
            # Plot and document mesh on the current refinement iteration
            postprocess.document_results(mesh, test_params, project_path, results_folder, results_info,
                                         ed_distr, md_distr, maj_distr, error, maj, i_eff_maj, min, i_eff_min,
                                         h_max_array, h_min_array)
            # Refine mesh
            #mesh = self.execute_refiniment_strategy(mesh, test_params, e_distr, m_distr)
            mesh = self.execute_refiniment_strategy(mesh, test_params, ed_distr, md_distr)

            # Update functional spaces, BC, and stiffness/mass matrices
            V, VV, V_exact, VV_exact, H_div = problem.functional_spaces(mesh,
                                                                           test_params,
                                                                           dim)
            # Update boundary faces
            boundary_facets, ds, dS = self.construct_boundary_faces(mesh)

            i = i + 1

        decay_result_folder = postprocess.create_decay_tag(test_params, project_path, results_folder, results_info)

        postprocess.document_errors_decay(float(dim), test_params, decay_result_folder,
                                          DOFs[0:i], h_min_array[0:i],
                                          e_H1_array[0:i], e_L2_array[0:i],
                                          maj_array[0:i], min_array[0:i], i_eff_maj_array[0:i],
                                          e_L_array[0:i], e_id_array[0:i])

    def execute_refiniment_strategy(self, mesh, test_params, e_distr, m_distr):
        # Define the refinement strategy
        if test_params['refinement_tag'] == "adaptive":
            # Define the distribution (error or majorant) upon which the refinement is based
            if test_params['refinement_criteria_tag'] == 'error':
                distr = e_distr
            elif test_params['refinement_criteria_tag'] == 'majorant':
                distr = m_distr
                # distr = maj_distr
            # elif test_params['refinement_criteria_tag'] == 'e-dwr':
            #    distr = e_dwr_distr
            # elif test_params['refinement_criteria_tag'] == 'residual':
            #    distr = residual_distr

            # Run the marking procedure
            theta = 0.2  # percentage of the marked cells

            if test_params['marking_tag'] == 'average':
                marking = estimates.averaged_marking(mesh, distr)
            elif test_params['marking_tag'] == 'predefined':
                marking = estimates.predefined_amount_of_elements_marking(mesh, distr, theta)
            elif test_params['marking_tag'] == 'bulk':
                marking = estimates.bulk_marking(mesh, distr, theta)

            # Refine mesh based on the marking criteria
            mesh = refine(mesh, marking)


        elif test_params['refinement_tag'] == "uniform":

            # Refinement for the comparison with theta < 1.0
            # cell_markers = CellFunction("bool", mesh)
            # cell_markers.set_all(True)
            # mesh = refine(mesh, cell_markers)

            mesh = refine(mesh)

        return mesh

    # Refinement for the optimal convergence

# Define refinement strategies
adaptive_tag = "adaptive"
uniform_tag = "uniform"

# Define marking strategies
average_tag = "average"
bulk_tag = "bulk"
predefined_tag = "predefined"

# Define the criteria of the refinement
error_tag = "error"
majorant_tag = "majorant"
e_dwr_tag = "e-dwr"
residual_tag = "residual"

# Define solving strategy
with_refinement_tag = "with-refinement"
without_refinement_tag = "without-refinement"

# Dictionary of tests
tests = {1: tests.quadratic_polynomial_solution_2d_t,
         2: tests.quadratic_polynomial_solution_1d_t,
         3: tests.quadratic_polynomial_solution_1d_t_2,
         4: tests.quadratic_polynomial_solution_2d_t_2,
         5: tests.sinusoidal_polynomial_solution_1d_t,
         6: tests.solution_with_singularities_2d_t,
         7: tests.sinusoidal_polynomial_solution_1d_t_on_long_rectangular_domain,
         8: tests.example_steinbach_paper_1d_t_c_1_k_1_cH_1,
         9: tests.example_steinbach_paper_1d_t_c_1_k_1_cH_10,
         10: tests.example_tutorial_example_3_a_x_cH_10,
         11: tests.example_steinbach_paper_1d_t_c_1_k_1_cH_1_with_u0,
         12: tests.example_tutorial_example_2,
         13: tests.example_steinbach_paper_1d_t_c_1_k_1_cH_100_with_u0,
         14: tests.example_steinbach_paper_1d_t_c_100,
         15: tests.example_tutorial_example_2_cH_10,
         16: tests.example_tutorial_example_2_cH_20,
         17: tests.example_tutorial_example_2_cH_1,
         18: tests.example_tutorial_example_2_cH_100,
         19: tests.quadratic_polynomial_solution_1d_t_2_A_heterogenous,
         20: tests.quadratic_polynomial_solution_1d_t_2_A_heterogenous_quadr,
         21: tests.sin_solution_1d_t_2_A_heterogenous_lin,
         22: tests.sin_solution_1d_t_2_A_heterogenous_quadr,
         23: tests.sin_solution_1d_t_2_A_heterogenous_sin,
         24: tests.fokker_plank_example_2_1d,
         25: tests.fokker_plank_example_3_sigma_1eminus1_1d,
         26: tests.linear_trigonometric_solution_1d_rho_1minus3_t,
         27: tests.quadratic_polynomial_solution_conv_react_sigma_1_1d_t,
         28: tests.example_tutorial_example_3_a_const_lambda_sigma_1eminus1_cH_2,
         29: tests.example_tutorial_example_3_a_const_lambda_sigma_5eminus2_cH_2,
         30: tests.example_tutorial_example_3_a_const_lambda_sigma_1_cH_2,
         31: tests.quadratic_polynomial_solution_1d_t_with_reference_solution_calculation,
         32: tests.quadratic_polynomial_solution_1d_t_with_exact_solution_calculation,
         33: tests.quadratic_polynomial_solution_1d_t_2_with_reference_solution,
         34: tests.example_tutorial_example_3_a_const_lambda_sigma_1_cH_2_with_reference,
         35: tests.example_tutorial_example_3_a_const_lambda_sigma_1eminus2_cH_1_with_reference,
         36: tests.example_tutorial_example_3_a_const_lambda_sigma_1eminus2_cH_1,
         37: tests.example_tutorial_example_4_a_linear_lambda_sigma_1eminus2_cH_1_Neumann_right,
         38: tests.example_tutorial_example_5_a_linear_lambda_sigma_1eminus2_cH_1_Neumann_pure,
         39: tests.example_tutorial_example_6_a_linear_lambda_sigma_1eminus2_cH_1_Neumann_pure_zero,
         40: tests.example_tutorial_example_7_a_linear_lambda_sigma_1eminus3_cH_1_Neumann_pure_zero,
         41: tests.example_tutorial_example_7_a_linear_lambda_sigma_1e3_cH_1_Neumann_pure_zero,
         42: tests.example_tutorial_example_8_a_const_cH_1_Neumann_pure_zero}
         #,
         #31: tests.example_tutorial_example_3_a_const_lambda_sigma_5eminus3_cH_2}

# Set the number of the test and call for the problem data
test_num = 40
u_expr, grad_u_expr, f_expr, \
A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps, \
u0_expr, uD_expr, uN_expr, T, domain, dim, \
solution_tag, material_tag, pde_tag = tests[test_num]()

# Init problem parameter
# s
test_params = {'refinement_tag': adaptive_tag,
                  'marking_tag': bulk_tag,
                  'percentage_value': 0.6,
                  'refinement_criteria_tag': error_tag,
                  'solving_strategy_tag': "with-refinement",
                  'number_of_refinements': 10,
                  'nx0': 8,
                  'nx1': 8,
                  'nx2': 8,
                  'nt':  16,
                  'res':  2,
                  'v_approx_order': 1,
                  'flux_approx_order': 2,
                  'v_exact_approx_order': expression_degree,
                  'expression_degree': expression_degree,
                  'solution_tag': solution_tag,
                  'material_tag': material_tag,
                  'pde_tag': pde_tag,
                  'majorant_optimization_iterations': 3,
                  'error_estimates': True,
                  'MAJORANT_OPTIMIZE': True,
                  'PLOT': True,
                  'STABILIZE': True,
                  'LOG_TO_FILE': False}

test = TestEstimates(test_num, u_expr, grad_u_expr, f_expr, 
                     A_expr, min_eig_A, lambda_expr, a_expr, c_H, eps,
                     u0_expr, uD_expr, uN_expr, T, domain, dim)
test.test_estimates(test_params)




