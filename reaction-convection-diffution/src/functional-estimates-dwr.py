__author__ = "Svetlana Matculevich <svetlana.v.matculevich@jyu.fi>"
__date__ = "2013-10-11"
__copyright__ = "Copyright (C) 2013 Svetlana Matculevich"
__license__ = ""

# Last changed: 2013-10-13

from dolfin import *
from dolfin.cpp.mesh import UnitSquareMesh, Mesh, UnitIntervalMesh, UnitCubeMesh, MeshFunction, cells, CellFunction, \
    Point
import time
import tests, problem, estimates, postprocess
from problem import Grad, Div, inverse
from mshr import *

def execute_refinement_strategy(problem_params, mesh, e_distr, md_distr,
                                project_path, results_folder, init_data_info):
    # error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_max_k, h_min_k, t_T, k, nt,

        # Define the refinement strategy
    if problem_params['refinement_tag'] == "adaptive":
        # Define the distribution (error or majorant) upon which the refinement is based
        if problem_params['refinement_criteria_tag'] == 'error':
            distr = e_distr
        elif problem_params['refinement_criteria_tag'] == 'majorant':
            distr = md_distr

        if problem_params['marking_tag'] == 'average':
            marking = estimates.averaged_marking(mesh, distr)
            refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking'

        elif problem_params['marking_tag'] == 'predefined':
            marking = estimates.predefined_amount_of_elements_marking(mesh, distr, problem_params['percentage_value'])
            refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking' + '-proc-%d' % (problem_params['percentage_value'] * 100)

        elif problem_params['marking_tag'] == 'bulk':
            marking = estimates.bulk_marking(mesh, distr, problem_params['percentage_value'])
            refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking' + '-bulk-%d' % (problem_params['percentage_value'] * 100)
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

#-------------------------------------------------------------------------------#
def solve_problem_nd_t(mesh, V, VV, V_exact, VV_exact, H_div,
                       f, A, invA, adjA, A_expr, lambda_1, lmbd, a,
                       u_e, uD, uD_boundary,
                       dim, domain, C_FD, nx0, nx1, nx2,
                       test_num,
                       results_folder,
                       problem_params):

    # Initialize the project path
    project_path = postprocess.get_project_path()

    if problem_params['solving_strategy_tag'] == "with-refinement":
        ref_num = 6
    else:
        ref_num = 1

    ref_num_array = postprocess.allocate_array(ref_num)
    h_max_array = postprocess.allocate_array(ref_num)
    h_min_array = postprocess.allocate_array(ref_num)


    e_array = postprocess.allocate_array(ref_num)

    maj_array = postprocess.allocate_array(ref_num)
    min_array = postprocess.allocate_array(ref_num)

    i_eff_maj_array = postprocess.allocate_array(ref_num)
    i_eff_min_array = postprocess.allocate_array(ref_num)

    for i in range(0, ref_num):

        print(" ")
        print("%------------------------------------------------------------------------------------------------------")
        print("% Solving the problem ")
        print("%------------------------------------------------------------------------------------------------------")
        print(" ")

        # Compute approximate solution
        u, delta = problem.solve_convection_reaction_diffusion(V, f, A, lambda_1, lmbd, problem.interpolate_vector_function(a, dim, V_exact, VV_exact), uD, uD_boundary, mesh, dim)
        if problem_params['solution_tag'] == "reference-solution":
            # Compute reference solution
            u_e, delta = problem.solve_convection_reaction_diffusion(V_exact, f, A, lambda_1, lmbd, problem.interpolate_vector_function(a, dim, V_exact, VV_exact), uD, uD_boundary, mesh, dim)

        # Plot approximate solution
        #if i == 0:
        #    p = plot(u)
        #else:
        #    p.plot(u)

        postprocess.plot_function(u, mesh, dim, project_path + results_folder + 'u-ref-%d' % i)
        
        # Calculate error
        error, var_grad_e, var_lmbd_diva_e = estimates.error_norm(u, u_e, lmbd, A, problem.interpolate_vector_function(a, dim, V_exact, VV_exact), V_exact, mesh, dim)

        # L2-projection of grad u to Hdiv space
        y = project(A * Grad(u, dim), H_div)
        #y = interpolate(A * Grad(u, dim), H_div)
        MAJORANT_OPTIMIZE = 1

        # Exact flux (for the reference)
        #y = project(A * Grad(u_e, dim), H_div)
        #MAJORANT_OPTIMIZE = 0

        # Calculate majorant
        maj, y, beta, md, mf, var_m_d, var_m_f_w_opt, \
        majorant_reconstruction_time, majorant_minimization_time = \
            estimates.majorant_nd(u, y, H_div, f, A, a, lambda_1, lmbd, error, mesh, dim, domain, C_FD, MAJORANT_OPTIMIZE)
        # Output the results of upper bound reconstructions
        estimates.output_result_error_and_majorant(error, maj, sqrt(maj / error))

        # Construct error and majorant distribution
        ed_distr, delta_e_distr, e_distr, md_distr, mf_distr, maj_distr, ed, md = \
            estimates.error_majorant_distribution_nd(mesh, dim, V_exact,
                                                     var_grad_e, var_lmbd_diva_e, var_m_d, var_m_f_w_opt, beta)

        # Calculate minorant
        #min, phi = estimates.minorant(u, delta, mesh, V_exact, u_0, u0_boundary, f, A, a, lmbd, dim, error)
        # Output the results of lower bound reconstructions
        #estimates.output_result_minorant(min, min/error)
        
        '''
        md_CG0_distr, mdf_CG0_distr, e_CGO_distr = estimates.majorant_distribution_DG0(mesh, f, lmbd, A, invA, u, e_form, y, beta, C_FD, dim) 

        residual_CG0_distr, e_DWR_CG0_distr, J_e_CG0_distr = estimates.get_indicators_CG0(mesh, V, V_exact, f, A, adjA, lmbd, u_0, u0_boundary, u, u_e,
                                                                                          e_form, sqrt(e_d), dim)

        # Plot histograms with obtained indicators
        postprocess.plot_histogram(mesh, J_e_CG0_distr, residual_CG0_distr,
                                   project_path + results_folder + 'je-residual-distr-hist' + results_info)
        postprocess.plot_histogram(mesh, J_e_CG0_distr, e_DWR_CG0_distr,
                                   project_path + results_folder + 'je-edwr-distr-hist' + results_info)
        # Plot histograms with obtained indicators from majorant
        postprocess.plot_histogram(mesh, J_e_CG0_distr, md_CG0_distr,
                                   project_path + results_folder + 'je-md-distr-hist' + results_info)
        postprocess.plot_histogram(mesh, J_e_CG0_distr, mdf_CG0_distr,
                                   project_path + results_folder + 'je-mdf-distr-hist' + results_info)
        '''

        # Update the arrays of data with respect to the refinement cycle
        ref_num_array[i] = i + 1
        e_array[i] = error
        maj_array[i] = maj
        i_eff_maj_array[i] = sqrt(maj / error)
        h_max_array[i] = mesh.hmax()
        h_min_array[i] = mesh.hmin()

        #min_array[i] = min
        #i_eff_min_array[i] = sqrt(min / error)

        num_cells = mesh.num_cells()
        num_vertices = mesh.num_vertices()

        # Construct the tag with problem information
        results_info = postprocess.construct_result_tag(test_num, i, nx0, nx1, nx2, num_cells, num_vertices)

        # Plot histogram with the error and majorant distribution
        postprocess.plot_histogram(mesh, ed_distr, md_distr,
                                   project_path + results_folder + 'e-maj-distr-hist' + results_info)

        # Refine mesh
        mesh, mat_file_tag = execute_refinement_strategy(problem_params, mesh, e_distr, md_distr,
                                                         project_path, results_folder, results_info)

        # Save the results in mat-file
        postprocess.save_results_to_mat_file(ed_distr, md_distr, e_array, maj_array, i_eff_maj_array,
                                             mat_file_tag)

        if problem_params['refinement_tag'] == "adaptive":
            # Plot result mesh
            if dim == 2:
                postprocess.plot_mesh(mesh, project_path + results_folder + 'mesh' + results_info)
                postprocess.plot_carpet_2d(mesh, ed, project_path + results_folder + 'carpet-error' + results_info)
                postprocess.plot_carpet_2d(mesh, md, project_path + results_folder + 'carpet-majorant' + results_info)

            elif dim == 3:
                postprocess.plot_mesh_3d(mesh, project_path + results_folder + 'initial-mesh')

        if i == ref_num - 1:

            # Define the refinement_tag
            if problem_params['refinement_tag'] == "adaptive":
                refinement_tag = problem_params['refinement_criteria_tag'] + '-' + \
                                 problem_params['marking_tag'] + '-marking'
            elif problem_params['refinement_tag'] == "uniform":
                refinement_tag = problem_params['refinement_tag']

            postprocess.plot_decay_error_majorant(problem_params["v_approx_order"], h_min_array, e_array, maj_array,
                                                  project_path + results_folder +
                                                  "error-maj-v-P%d-" % problem_params["v_approx_order"] +
                                                  refinement_tag + results_info)

            postprocess.save_decay_to_mat_file(h_min_array, e_array, maj_array,
                                               project_path + results_folder +
                                               "maj-decay-v-P%d" % problem_params["v_approx_order"])
            postprocess.plot_decay_majorant_v_of_deg(h_min_array, maj_array, problem_params["v_approx_order"],
                                                     project_path + results_folder +
                                                     "maj-decay-v-P%d" % problem_params["v_approx_order"] +
                                                     results_info)

        if ref_num > 1 and i < ref_num - 1:

            # Update functional spaces, BC, and stiffness/mass matrices
            V, VV, V_exact, VV_exact, H_div = problem.functional_spaces(mesh,
                                                                        problem_params["v_approx_order"],
                                                                        problem_params["flux_approx_order"],
                                                                        dim)
            if problem_params['material_tag'] == "material-changing": # over the domain
                # Define A if it is changing
                A = problem.construct_from_mesh_functions(dim, A_expr, mesh)
        '''
        #tag = 'E-DWR-'
        #tag = 'eta-'
        #tag = 'error-'
        #tag = 'm-d-'
        tag = 'm-df-'


        eta_distr, E_DWR_distr, J_e_distr, m_d_distr, m_df_distr  = \
            compare_error_indicators(mesh, V, V_star, V_e, f, u_0, u0_boundary, u,
                                     interpolate(u_e, V_e), interpolate(grad_u_e, VV_e),
                                     y, beta, test_num, tag)



        cells_num = mesh.num_cells()
        N_array[i] = cells_num
        e_array[i] = sum(J_e_distr)

        if tag == 'eta-':
            distr = eta_distr
        elif tag == 'E-DWR-':
            distr = E_DWR_distr
        elif tag == 'error-':
            distr = J_e_distr
        elif tag == 'm-d-':
            distr = m_d_distr
        elif tag == 'm-df-':
            distr = m_df_distr
            '''
#----------------------------------------------------------------------------------------------------------------------#
# Class TestEstimates
#----------------------------------------------------------------------------------------------------------------------#
class TestEstimates():

    def __init__(self, test_num, u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim):

        # Initialize the problem data
        self.u_expr = u_expr
        self.grad_u_expr = grad_u_expr
        self.lambda_expr = lambda_expr
        self.f_expr = f_expr
        self.A_expr = A_expr
        self.lambda_1 = lambda_1
        self.a_expr = a_expr
        self.u0_boundary = u0_boundary
        self.uD_expr = uD_expr
        self.dim = dim
        self.domain = domain
        self.test_num = test_num

    def convert_problem_data_to_expressions(self, mesh, problem_params):

        if problem_params['material_tag'] == "material-changing": # over the domain
            # Define A if it is changing
            A = problem.construct_from_mesh_functions(dim, A_expr, mesh)

        elif problem_params['material_tag'] == "material-constant":
            # Convert matrix to UFL form
            if dim == 1:
                A = Expression(A_expr[0])
            else:
                A = as_matrix(A_expr)
        # Inverse and agjoint matrices
        invA = inverse(A, self.dim)
        adjA = det(A) * inverse(A, self.dim)

        if dim == 1:
            a = Expression(a_expr[0])
        elif dim >= 2:
            if dim == 2:
                a = Expression((a_expr[0], a_expr[1]))
            elif dim == 3:
                a = Expression((a_expr[0], a_expr[1], a_expr[2]))
   
        # Define exact solution expression
        if problem_params['solution_tag'] == "reference-solution": # over the domain
            u_e = "0"
            grad_u_e = []
        elif problem_params['solution_tag'] == "predefined-solution": # over the domain
            u_e = Expression(self.u_expr)
            grad_u_e = Expression(self.grad_u_expr)

        # Define Dirichlet BC function
        uD = Expression(self.uD_expr)

        # Define input data right-hand side
        f = Expression(self.f_expr)
        lmbd = Expression(self.lambda_expr)

        return f, A, invA, adjA, lmbd, a, uD, u_e, grad_u_e


    def test_estimates(self, problem_params):

        # Define the project path
        project_path = postprocess.get_project_path()

        # Define the mesh based on the domain geometry
        #--------------------------------------------------------------------------------------------------------------#
        if self.domain == "unit-domain":
            nx0 = 10
            nx1 = 10
            nx2 = 10

            if self.dim == 1:
                # Create mesh on the unit interval
               mesh = UnitIntervalMesh(nx0)
            elif self.dim == 2:
                # Create mesh on the unit square
                mesh = UnitSquareMesh(nx0, nx1)
            elif self.dim == 3:
                # Create mesh on the unit cube
                mesh = UnitCubeMesh(nx0, nx1, nx2)
        elif self.domain == "circle-domain":
            nx0 = 0
            nx1 = 0
            nx2 = 0

            res = 4
            circle = Circle(Point(0.0, 0.0), 2.0)
            mesh = generate_mesh(circle, res)
            #plot(mesh)

        elif self.domain == "l-shape-domain":
            nx0 = 0
            nx1 = 0
            nx2 = 0

            # Define the data path
            file_path = '/data/%dd/test-%d/' % (self.dim, test_num)

            # 2d l-shaped domain mesh
            file_name = 'l-shape-diam-2sqrt2-33-vertices.xml'
            # 3d l-shaped domain mesh
            #file_name = 'l-shape-3d-twice-refined.xml'
            #file_name = 'l-shape-3d-once-refined.xml'
            #file_name = 'l-shape-3d-three-refined.xml'

            # Load mesh from the xml-file
            mesh = Mesh(project_path + file_path + file_name)

            # Plot and save the initial mesh
            if self.dim == 2:
                postprocess.plot_mesh(mesh, project_path + file_path + 'initial-mesh')
            elif self.dim == 3:
                postprocess.plot_mesh_3d(mesh, project_path + file_path + 'initial-mesh')
        elif self.domain == "ring-domain":

            nx0 = 0
            nx1 = 0
            nx2 = 0

            res = 4
            big_circle = Circle(Point(0.0, 0.0), 1.0)
            small_circle = Circle(Point(0.0, 0.0), 0.5)
            mesh = generate_mesh(big_circle - small_circle, res)
            plot(mesh, interactive=True)

        # Calculate the estimate for Freidrichs constant based on the domain
        C_FD = problem.calculate_CF_of_domain(domain, self.dim)


        # Check if exist and create folder for the results
        #--------------------------------------------------------------------------------------------------------------#
        results_folder_name = postprocess.construct_result_folder_name(self.dim, test_num, problem_params)
        postprocess.create_results_folder(results_folder_name)

        # Define mesh and functional spaces based on the mesh
        #--------------------------------------------------------------------------------------------------------------#
        V, VV, V_exact, VV_exact, H_div = \
            problem.functional_spaces(mesh,
                                      problem_params["v_approx_order"],
                                      problem_params["flux_approx_order"],
                                      self.dim)

        # Define problem data functions
        #--------------------------------------------------------------------------------------------------------------#
        f, A, invA, adjA, lmbd, a, uD, u_e, grad_u_e = self.convert_problem_data_to_expressions(mesh, problem_params)

        # Output to console the problem data
        problem.output_problem_characteristics(self.test_num, self.u_expr, self.grad_u_expr, self.f_expr, self.A_expr, self.lambda_expr, self.a_expr, self.uD_expr,
                                               self.domain, self.dim, mesh.num_cells(), mesh.num_vertices(),
                                               problem_params["v_approx_order"], problem_params["flux_approx_order"])

        t0_problem = time.clock()
        solve_problem_nd_t(mesh, V, VV, V_exact, VV_exact, H_div,
                           f, A, invA, adjA, self.A_expr, lambda_1, lmbd, a,
                           u_e,
                           uD, self.u0_boundary,
                           self.dim, self.domain, C_FD,
                           nx0, nx1, nx2,
                           self.test_num,
                           results_folder_name,
                           problem_params)
        t1_problem = time.clock()
        print("time = %d" %(t1_problem - t0_problem))

# Define domain type
unit_domain_tag = "unit-domain"
l_shape_domain_tag = "l-shape-domain"
pi_shape_domain_tag = "pi-shape-domain"
circle_domain_tag = "circle-domain"

# Define solving strategy
with_refinement_tag = "with-refinement"
without_refinement_tag = "without-refinement"

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
material_changing_tag = "material-changing"
material_constant_tag = "material-constant"

# Define whether solution is known or not
predefined_solution_tag = "predefined-solution"
reference_solution_tag = "reference-solution"

# Functional spaces parameters
v_deg = 1
y_deg = 2

# Dictionary of tests
tests = {1: tests.quadratic_polynomial_solution_1d_sigma_1, # checked
         2: tests.quadratic_polynomial_solution_1d_sigma_1000, # checked
         4: tests.linear_trigonometric_solution_1d_rho_1minus3, # checked
         5: tests.linear_trigonometric_solution_1d_rho_1, #checked
         6: tests.linear_trigonometric_solution_1d_rho_1000, # checked
         7: tests.quadratic_polynomial_solution_2d_sigma_x_1_sigma_y_5, # checked
         8: tests.quadratic_polynomial_solution_2d_sigma_x_1minus2_sigma_y_5, # checked
         9: tests.linear_trigonometric_solution_2d_rho_1000, # checked
         10: tests.quadratic_polynomial_solution_1d_sigma_1minus2,# checked
         11: tests.quadratic_polynomial_solution_1d_lmdb_zero, # checked
         12: tests.quadratic_polynomial_solution_2d_with_A, # checked
         13: tests.quadratic_polynomial_solution_2d_with_A_sigma_5minus2, # checked
         14: tests.quadratic_polynomial_solution_3d_with_A, # checked, 3 iterations
         15: tests.quadratic_polynomial_solution_1d_sigma_1_epsilon_1eminus2, # checked
         16: tests.quadratic_polynomial_solution_1d_lmdb_zero_eps_1minus2, # checked
         17: tests.quadratic_polynomial_solution_1d_lmdb_zero_a_linear, # checked
         18: tests.quadratic_polynomial_solution_1d_lmdb_zero_a_cubic, # checked
         19: tests.quadratic_polynomial_solution_1d_lmdb_zero_a_cubic_epsilon_1minus2, # checked
         20: tests.classic_example_1d_b_1_eps_1, # checked
         21: tests.classic_example_1d_b_1_eps_1eminis2, # checked
         22: tests.f_1_uD_polynomial_1d_sigma_1000, # checked
         23: tests.f_1_uD_polynomial_1d_sigma_1, # checked
         24: tests.f_1_uD_polynomial_1d_sigma_1minus2, # checked
         25: tests.f_1_uD_polynomial_1d_lmdb_zero, # checked
         26: tests.f_trigonometric_uD_0_1d_rho_1minus3, # checked
         27: tests.f_trigonometric_uD_0_1d_rho_1000, # checked
         28: tests.f_trigonometric_uD_0_1d_rho_1, # checked
         29: tests.f_trigonometric_uD_0_2d_rho_1000, # checked
         30: tests.f_1_uD_0_2d_with_A,
         31: tests.f_1_uD_0_2d_with_A_sigma_x_5minus2_sigma_y_5,
         32: tests.f_1_uD_0_2d_with_changing_A_lmdb_0,
         33: tests.solution_with_singularities_2d,
         34: tests.quadratic_polynomial_solution_3d_with_A,
         35: tests.f_1_uD_zero_1d_lmdb_zero_a_const,
         36: tests.f_1_uD_polynomial_1d_lmdb_zero_a_const_eps_1minus3,
         37: tests.f_1_uD_zero_1d_lmdb_zero_a_cubic_eps_1minus3,
         38: tests.f_1_uD_zero_2d_lmdb_zero_a_const_eps_1minus2,
         39: tests.f_1_uD_zero_2d_lmdb_zero_a_const_eps_1minus3,
         40: tests.classic_example_1d_b_3_eps_1eminis2,
         41: tests.example_3_book_stynes_2d,
         42: tests.example_1_book_stynes_2d,
         43: tests.example_2_eps_1eminus2book_stynes_2d,
         44: tests.example_2_eps_1eminus3book_stynes_2d}

# Set the number of the test and call for the problem data
test_num = 44
u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim, solution_tag, material_tag = tests[test_num]()


# Init problem parameters
problem_params = {'refinement_tag': uniform_tag,
                  'marking_tag': average_tag,
                  'percentage_value': 0.4,
                  'refinement_criteria_tag': majorant_tag,
                  'solving_strategy_tag': with_refinement_tag,
                  'v_approx_order': v_deg,
                  'flux_approx_order': y_deg,
                  "solution_tag": solution_tag,
                  "material_tag":material_tag}

# Define Dirichlet boundary condition
def u0_boundary(x, on_boundary):
    return on_boundary

test = TestEstimates(test_num, u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, domain, dim)
test.test_estimates(problem_params)




