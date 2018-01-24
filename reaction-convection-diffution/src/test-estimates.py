__author__ = "Svetlana Matculevich <svetlana.v.matculevich@jyu.fi>"
__date__ = "2013-10-11"
__copyright__ = "Copyright (C) 2013 Svetlana Matculevich"
__license__ = ""

# Last changed: 2013-10-13

from dolfin import *
from mshr import *
from dolfin.cpp.mesh import UnitSquareMesh, Mesh, UnitIntervalMesh, UnitCubeMesh, CellFunction, Point, SubDomain, \
    FacetFunction
import time
import tests, problem, estimates, postprocess
from problem import Grad, inverse
from dolfin.cpp.common import set_log_active, set_log_level
from dolfin.cpp.function import near


CRITICAL  = 50 # errors that may lead to data corruption and suchlike
ERROR     = 40 # things that go boom
WARNING   = 30 # things that may go boom later
INFO      = 20 # information of general interest
PROGRESS  = 16 # what's happening (broadly)
TRACE     = 13 # what's happening (in detail)
DBG       = 10  # sundry
set_log_level(20)

parameters["plotting_backend"] = "matplotlib"
#parameters["refinement_algorithm"] = "plaza_with_different_faces"
#parameters["refinement_algorithm"] = "plaza_with_parent_faces"
#parameters["form_compiler"]["quadrature_degree"] = 14
#parameters["form_compiler"]["optimize"] = True
#parameters["form_compiler"]["representation"] = 'quadrature'

expression_degree = 4


# Define on boundary function
def boundary(x, on_boundary):
    return on_boundary

# Class of the Dirichlet BC
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary

class LeftBoundary(AutoSubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0)

class BottomBoundary(AutoSubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[1], 0)

#----------------------------------------------------------------------------------------------------------------------#
# Class TestEstimates
#----------------------------------------------------------------------------------------------------------------------#
class TestEstimates():

    def __init__(self, test_num,
                       u_expr, grad_u_expr, f_expr, 
                       A_expr, lambda_1, lambda_expr, a_expr, uD_expr, uN_expr,
                       domain, dim):

        # Initialize the problem data
        self.u_expr = u_expr
        self.grad_u_expr = grad_u_expr
        self.lambda_expr = lambda_expr
        self.f_expr = f_expr
        self.A_expr = A_expr
        self.lambda_1 = lambda_1
        self.a_expr = a_expr
        self.boundary = boundary
        self.uD_expr = uD_expr
        self.uN_expr = uN_expr
        self.dim = dim
        self.domain = domain
        self.test_num = test_num

    def calculate_trace_constant(self):

        # No Neumann part of the boundary
        C_Ntr = 0.0
        # Mark different parts of the domain for non-homogeneous Dirichlet BC
        if test_num == 48:
            # Caclulated in the matlab
            C_Ntr = 1.47

        return C_Ntr

    def construct_boundary_faces(self, mesh):

        # Get boundary facets function
        facet_funcs = FacetFunction('size_t', mesh)
        # Dirichlet part of the boudary is marked by 0
        dirichlet_marker = 0
        facet_funcs.set_all(dirichlet_marker)

        neumann_marker = dirichlet_marker + 1

        # Mark different parts of the domain for non-homogeneous Dirichlet BC
        if test_num == 45:
            # Create parts of the boundary as subdomains
            bottom = AutoSubDomain(lambda x: near(x[1], 0))
            # Mark part of the domain with non-zero Dirichlet condition
            bottom.mark(facet_funcs, dirichlet_marker + 1)

        elif test_num == 48:
            # Create parts of the boundary as subdomains
            right = AutoSubDomain(lambda x: near(x[0], 1))
            # Mark part of the domain with non-zero Dirichlet condition
            right.mark(facet_funcs, neumann_marker)

        elif test_num == 70:
            # Create parts of the boundary as subdomains
            right = AutoSubDomain(lambda x: near(x[0], 1))
            # Mark part of the domain with non-zero Dirichlet condition
            right.mark(facet_funcs, neumann_marker)

            # Create parts of the boundary as subdomains
            left = AutoSubDomain(lambda x: near(x[0], 0.0))
            # Mark part of the domain with non-zero Dirichlet condition
            left.mark(facet_funcs, neumann_marker)

            # Create parts of the boundary as subdomains
            top = AutoSubDomain(lambda x: near(x[1], 1.0))
            # Mark part of the domain with non-zero Dirichlet condition
            top.mark(facet_funcs, neumann_marker)

            # Create parts of the boundary as subdomains
            bottom = AutoSubDomain(lambda x: near(x[1], 0.0))
            # Mark part of the domain with non-zero Dirichlet condition
            bottom.mark(facet_funcs, neumann_marker)

        ds = Measure('ds')[facet_funcs] # outer facets
        dS = Measure('dS')[facet_funcs] # inner facets

        return facet_funcs, ds, dS

    # Convert given data expressions form string format to the UFL Expressions
    def convert_problem_data_to_expressions(self, mesh, problem_params):
        """
        :param mesh: mesh-discretization of the domain
        :param dim: dimension of the problem
        :return: given data as functions
                 f, A, invA, adjA, lmbd, a, uD, u_e, grad_u_e
        """

        # Define diffusion matrix A if it is changing on different parts of the domain
        if problem_params['material_tag'] == "material-changing":
            # Construct A from the mesh function
            A = problem.construct_from_mesh_functions(dim, A_expr, mesh)

        # Define diffusion matrix A if it is constant over the whole domain
        elif problem_params['material_tag'] == "material-constant":
            # Convert matrix to UFL form depanding on
            if dim == 1:
                # Convert string expression to UFLExpression
                A = Expression(A_expr[0], degree=expression_degree)
            else:
                # Convert matrix to UFLMatrix
                A = as_matrix(A_expr)

        # Construct inverse and adjoint matrices for the diffusion matrix A
        invA = inverse(A, A_expr, self.dim)
        adjA = det(A) * invA

        # Construct expression for the convection/advection vector field
        if dim == 1:
            a = Expression(a_expr[0], degree=expression_degree + 2)
        elif dim >= 2:
            if dim == 2:
                a = Expression((a_expr[0], a_expr[1]), degree=expression_degree)
            elif dim == 3:
                a = Expression((a_expr[0], a_expr[1], a_expr[2]), degree=expression_degree)

        # Construct expression for the exact solution expression
        if problem_params['solution_tag'] == "reference-solution":
            u_e = "0"
            grad_u_e = []
        elif problem_params['solution_tag'] == "predefined-solution":
            u_e = Expression(self.u_expr, degree=expression_degree)
            grad_u_e = Expression(self.grad_u_expr, degree=expression_degree)

        # Construct the Dirichlet BC depending on the concrete example
        if test_num == 45:
            # Define Dirichlet BC out of two functions
            uD = [Expression(self.uD_expr[0], degree=expression_degree), \
                  Expression(self.uD_expr[1], degree=expression_degree)]

        else:
            # Define Dirichlet BC function
            uD = Expression(self.uD_expr, degree=expression_degree)

        uN = Expression(self.uN_expr, degree=problem_params["v_approx_order"])

        # Construct the expression for the right-hand side
        f = Expression(self.f_expr, degree=expression_degree)

        # Construct the expression for the reaction function
        lmbd = Expression(self.lambda_expr, degree=expression_degree)

        return f, A, invA, adjA, lmbd, a, uD, uN, u_e, grad_u_e

    # Function to generate the mesh
    def construct_mesh(self, nx0, nx1, nx2, res, project_path):
        """
        :param nx0, nx1, nx2: mesh-sizes wrt to X, Y, Z directions
        :param res: mesh resolution if defined with mesh_generator function
        :param project_path: path of the project to load the mesh
        :return: mesh:discretization of the domain
        """
        # --------------------------------------------------------------------------------------------------------------#
        if self.domain == "unit-domain":

            if self.dim == 1:
                # Create mesh on the unit interval
                mesh = UnitIntervalMesh(nx0)

            elif self.dim == 2:
                # Create mesh on the unit square
                mesh = UnitSquareMesh(nx0, nx1)

            elif self.dim == 3:
                # Create mesh on the unit cube
                mesh = UnitCubeMesh(nx0, nx1, nx2)
        # --------------------------------------------------------------------------------------------------------------#
        elif self.domain == "circle-domain":
            # Define the 'resolution' of the mesh
            res = 4

            # Define the radius of the circle
            rad = 2.0

            # Create the circle geometry
            circle = Circle(Point(0.0, 0.0), rad)

            # Generate the mesh based on the geometry and the resolution of the mesh
            mesh = generate_mesh(circle, res)
        # --------------------------------------------------------------------------------------------------------------#
        elif self.domain == "l-shape-domain":
            # Define the data path depending on the problem parameters
            file_path = '/data/%dd/test-%d/' % (self.dim, test_num)

            # 2d l-shaped domain mesh
            file_name = 'l-shape-diam-2sqrt2-33-vertices.xml'

            # 3d l-shaped domain mesh
            # file_name = 'l-shape-3d-twice-refined.xml'
            # file_name = 'l-shape-3d-once-refined.xml'
            # file_name = 'l-shape-3d-three-refined.xml'

            # Load mesh from the xml-file
            mesh = Mesh(project_path + file_path + file_name)
        # --------------------------------------------------------------------------------------------------------------#
        elif self.domain == "ring-domain":

            # Define the 'resolution' of the mesh
            res = 4

            # Define the radius of the circle
            big_rad = 1.0
            small_rad = 0.5

            # Create the circle geometry
            big_circle = Circle(Point(0.0, 0.0), big_rad)
            small_circle = Circle(Point(0.0, 0.0), small_rad)

            # Generate the mesh based on the geometry and the resolution of the mesh
            mesh = generate_mesh(big_circle - small_circle, res)

        return mesh

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

    # Main function of the clas for testing the estimates
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
        mesh = self.construct_mesh(test_params["nx0"], test_params["nx1"], test_params["nx2"], test_params["res"], project_path)

        if test_params["PLOT"] == True:
            # Plot and save the initial mesh
            if self.dim == 2:   postprocess.plot_mesh(mesh, project_path + results_folder_name + 'initial-mesh')
            elif self.dim == 3: postprocess.plot_mesh_3d(mesh, project_path + results_folder_name + 'initial-mesh')

        boundary_facets, ds, dS = self.construct_boundary_faces(mesh)

        # Calculate the estimate for Freidrichs constant based on the domain
        C_FD = problem.calculate_CF_of_domain(domain, self.dim)
        C_Ntr = self.calculate_trace_constant()

        # Define functional spaces based on the mesh
        #--------------------------------------------------------------------------------------------------------------#
        V, VV, V_exact, VV_exact, H_div = \
            problem.functional_spaces(mesh, test_params, self.dim)

        # Define problem data functions
        #--------------------------------------------------------------------------------------------------------------#
        f, A, invA, adjA, lmbd, a, uD, uN, u_e, grad_u_e = self.convert_problem_data_to_expressions(mesh, test_params)

        # Output to console the problem data
        problem.output_problem_characteristics(self.test_num, self.u_expr, self.grad_u_expr, self.f_expr, self.A_expr, self.lambda_expr, self.a_expr, self.uD_expr,
                                               self.domain, self.dim, mesh.num_cells(), mesh.num_vertices(),
                                               test_params["v_approx_order"], test_params["flux_approx_order"])

        # Run the solver with a posteriori error control
        # --------------------------------------------------------------------------------------------------------------#
        t0_problem = time.clock()
        self.solve_problem_nd(mesh, boundary_facets, ds, dS,
                                V, V_exact, VV_exact, H_div,
                                u_e, f, A, invA, adjA, lambda_1, lmbd, a,
                                uD, uN,
                                C_FD, C_Ntr,
                                project_path, results_folder_name,
                                test_params)

        t1_problem = time.clock()
        print("time = %d" %(t1_problem - t0_problem))

    # -------------------------------------------------------------------------------#
    def solve_problem_nd(self, mesh, boundary_facets, ds, dS,
                           V, V_exact, VV_exact, H_div,
                           u_e, f, A, invA, adjA, lambda_1, lmbd, a,
                           uD, uN,
                           C_FD, C_Ntr,
                           project_path, results_folder,
                           test_params):

        # Define the number of refinements
        ref_num = test_params['number_of_refinements']

        # Define arrays with data collected on the refinement steps
        ref_num_array = postprocess.allocate_array(ref_num + 1)
        h_max_array = postprocess.allocate_array(ref_num + 1)
        h_min_array = postprocess.allocate_array(ref_num + 1)
        dofs_array = postprocess.allocate_array(ref_num + 1)
        numcells_array = postprocess.allocate_array(ref_num + 1)
        numverts_array = postprocess.allocate_array(ref_num + 1)
        e_array = postprocess.allocate_array(ref_num + 1)
        e_l2_array = postprocess.allocate_array(ref_num + 1)
        e_linf_array = postprocess.allocate_array(ref_num + 1)

        e_min_array = postprocess.allocate_array(ref_num + 1)
        maj_array = postprocess.allocate_array(ref_num + 1)
        min_array = postprocess.allocate_array(ref_num + 1)
        i_eff_maj_array = postprocess.allocate_array(ref_num + 1)
        i_eff_min_array = postprocess.allocate_array(ref_num + 1)


        # Initialize the variables before the loop
        maj = 10e8
        i = 0


        # TO DO: Calculate the reference solution on the refined mesh
        '''
        if test_params['solution_tag'] == "reference-solution":
            mesh_ref = self.construct_mesh(2 ** ref_num * test_params["nx0"],
                                           2 ** ref_num * test_params["nx1"],
                                           2 ** ref_num * test_params["nx2"],
                                           2 ** ref_num * test_params["res"], project_path)
            # Update functional spaces, BC, and stiffness/mass matrices
            V_ref, VV_ref, V_ref_exact, VV_ref_exact, H_ref_div = problem.functional_spaces(mesh_ref,
                                                                            test_params["v_approx_order"],
                                                                            test_params["flux_approx_order"],
                                                                            dim)
            boundary_facets_ref, ds = self.construct_boundary_faces(mesh_ref)
            u_e, delta = problem.solve_convection_reaction_diffusion(V_ref_exact, f, A, lambda_1, lmbd,
                                                                     problem.interpolate_vector_function(a, dim, V_ref_exact, VV_ref_exact),
                                                                     self.boundary, uD, uN,
                                                                     mesh_ref, boundary_facets_ref,
                                                                     dim, test_num)
        '''

        while i <= ref_num: #and maj > test_params['accuracy_level']:

            print(" ")
            print("%-----------------------------------------------------------------------------------------------------------%")
            print " Refinement cycle # %d : DOFs = %d" %(i, len(V.dofmap().dofs()))
            print("%-----------------------------------------------------------------------------------------------------------%")
            print(" ")

            # Compute approximate solution u and stabilization parameter (in case of the convection dominated problem)

            u, delta = problem.solve_convection_reaction_diffusion(V, f, A, lambda_1, lmbd,
                                                                   problem.interpolate_vector_function(a, dim, V_exact,
                                                                                                       VV_exact),
                                                                   self.boundary, uD, uN,
                                                                   mesh, boundary_facets,
                                                                   self.dim,
                                                                   self.test_num, test_params)

            # In case of explicitely unknown exact solution, calculate it with high order order polynomials
            if test_params['solution_tag'] == "reference-solution":
                # Compute reference solution with high order polynomials
                u_e, delta = problem.solve_convection_reaction_diffusion(V_exact, f, A, lambda_1, lmbd,
                                                                         problem.interpolate_vector_function(a, dim,
                                                                                                             V_exact,
                                                                                                             VV_exact),
                                                                         self.boundary, uD, uN,
                                                                         mesh, boundary_facets,
                                                                         self.dim,
                                                                         self.test_num, test_params)

            if test_params["PLOT"] == True:
                #plot(u, interactive=True)
                # Plot approximate solution
                postprocess.plot_function(u_e, mesh, dim, project_path + results_folder + 'u-%d' % i)

            # Calculate error
            grad_e, e_l2, e_linf, delta_e, lambda_e, a_e, var_grad_e, var_lmbd_diva_e, var_lmbd_e, var_a_e = \
                estimates.error_norm(u, u_e, lmbd, A, invA,
                                     a,
                                     problem.interpolate_vector_function(a, dim, V_exact, VV_exact),
                                     V, V_exact, mesh, dim)
            # Define the error for the majorant
            error = grad_e + delta_e


            # Define the error for the minorant
            if test_params["pde_tag"] == "reaction-diffusion-pde" or test_params["pde_tag"] == "diffusion-pde":
                error_min = grad_e + delta_e
            elif test_params["pde_tag"] == "reaction-convection-diffusion-pde" or test_params["pde_tag"] == "convection-diffusion-pde":
                error_min = grad_e + lambda_e + a_e

            if test_params["error_estimates"] == True:

                # L2-projection of grad u to Hdiv space
                y = project(A * Grad(u, dim), H_div)

                # Test for the correctness of the code and the majorant for the y = A dot \grad u
                # y = interpolate(A * Grad(u, dim), H_div)
                #y = project(A * Grad(u_e, dim), H_div)

                # Contruct the majorant
                maj, y, beta, md, mf, var_m_d, var_m_f_w_opt, \
                majorant_reconstruction_time, majorant_minimization_time = \
                    estimates.majorant_nd(u, y, H_div, V_exact, VV_exact, f, A, invA, lambda_1, a, lmbd, error, mesh, dim, C_FD, C_Ntr, test_params)

                # Output the results of upper bound reconstructions
                estimates.output_result_error_and_majorant(sqrt(error), sqrt(maj), sqrt(maj / error))

                # Construct error and majorant distribution
                ed_distr, delta_e_distr, e_distr, md_distr, mf_distr, maj_distr, ed, md = \
                    estimates.error_majorant_distribution_nd(mesh, dim, V_exact,
                                                             var_grad_e, var_lmbd_diva_e, var_m_d, var_m_f_w_opt, beta)

                # Calculate minorant
                min, phi = estimates.minorant(u_e, u, mesh, boundary_facets, ds, V, V_exact, uD, uN, self.boundary, f, A, lmbd, a, dim, test_params)
                # Output the results of lower bound reconstructions
                if min < 0: min = 0
                estimates.output_result_minorant(sqrt(error_min), sqrt(min), sqrt(min/error_min))

                maj_array[i] = sqrt(maj)
                i_eff_maj_array[i] = sqrt(maj / error)

                min_array[i] = sqrt(min)
                i_eff_min_array[i] = sqrt(min / error_min)

                if test_params["pde_tag"] == "reaction-convection-diffusion-pde" \
                            or test_params["pde_tag"] == "convection-diffusion-pde":
                    e_min_array[i] = sqrt(error_min)
                else:
                    e_min_array[i] = e_array[i]

                '''
                #tag = 'E-DWR-'
                #tag = 'eta-'
                #tag = 'error-'
                #tag = 'm-d-'
                #tag = 'm-df-'

                eta_distr, E_DWR_distr, J_e_distr, m_d_distr, m_df_distr  = \
                    compare_error_indicators(mesh, V, V_star, V_e, f, u_0, boundary, u,
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

                '''
                md_CG0_distr, mdf_CG0_distr, e_CGO_distr = estimates.majorant_distribution_DG0(mesh, f, lmbd, A, invA, u, e_form, y, beta, C_FD, dim)

                residual_CG0_distr, e_DWR_CG0_distr, J_e_CG0_distr = estimates.get_indicators_CG0(mesh, V, V_exact, f, A, adjA, lmbd, u_0, boundary, u, u_e,
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
            e_array[i] = sqrt(error)
            e_l2_array[i] = e_l2
            e_linf_array[i] = e_linf
            h_max_array[i] = mesh.hmax()
            h_min_array[i] = mesh.hmin()

            num_cells = mesh.num_cells()
            num_vertices = mesh.num_vertices()
            num_dofs = len(V.dofmap().dofs())

            numcells_array[i] = num_cells
            numverts_array[i] = num_vertices
            dofs_array[i] = num_dofs

            # Construct the tag with problem information
            results_info = postprocess.construct_result_tag(test_num, i,
                                                            test_params["nx0"], test_params["nx1"], test_params["nx2"],
                                                            num_cells, num_vertices)

            # Plot histogram with the error and majorant distribution
            # postprocess.plot_histogram(mesh, ed_distr, md_distr,
            #                           project_path + results_folder + 'e-maj-distr-hist' + results_info)

            # For the refinement accept the last one change the mesh-dependent values like spaces and mesh functions
            if ref_num > 0 and i + 1 <= ref_num:

                # If the error estimates are constructed, then the mesh refinement can vary: uniform or adaptive
                if test_params["error_estimates"] == True:
                    # Refine mesh
                    mesh, mat_file_tag = self.execute_refinement_strategy(test_params, mesh, e_distr, md_distr,
                                                                          project_path, results_folder, results_info)
                    # Save the results in mat-file
                    postprocess.save_results_to_mat_file(ed_distr, md_distr, e_array, maj_array, i_eff_maj_array,
                                                         e_min_array, min_array,
                                                         mat_file_tag)

                    if i + 1 == ref_num: postprocess.save_mesh_to_xml_file(mesh,
                                                                           project_path + results_folder, results_info)

                    # If the refiniment strategy is adaptive, that plot the changes in the mesh, majorant, and the error
                    if test_params['refinement_tag'] == "adaptive":

                        if test_params["PLOT"] == True:
                            # Plot result mesh
                            if dim == 2:
                                postprocess.plot_mesh(mesh, project_path + results_folder + 'mesh' + results_info)

                                # Plot 'carpets with colored elements' dependent on the error and the majorant
                                postprocess.plot_carpet_2d(mesh, ed,
                                                           project_path + results_folder + 'carpet-error' + results_info)
                                postprocess.plot_carpet_2d(mesh, md,
                                                           project_path + results_folder + 'carpet-majorant' + results_info)
                            elif dim == 3:
                                postprocess.plot_mesh_3d(mesh, project_path + results_folder + 'initial-mesh')
                # If the error estimates aren't constructed, the mesh refinement is uniform
                else:
                    mesh = refine(mesh)

                # Update functional spaces, BC, and stiffness/mass matrices
                V, VV, V_exact, VV_exact, H_div = problem.functional_spaces(mesh,
                                                                            test_params, dim)
                if test_params['material_tag'] == "material-changing":  # over the domain
                    # Define A if it is changing
                    A = problem.construct_from_mesh_functions(dim, self.A_expr, mesh)

                boundary_facets, ds, dS = self.construct_boundary_faces(mesh)

                # Define the value of the maj for the loop criteria
                if test_params["error_estimates"] == True:
                    maj = maj_array[i]
                else:
                    maj = e_array[i]

            # Update the refinement number
            i += 1

        # Output the results
        decay_result_folder = postprocess.create_decay_tag(test_params, project_path, results_folder, results_info)

        postprocess.document_errors_decay(float(self.dim), test_params, decay_result_folder,
                                          dofs_array[0:i], h_min_array[0:i],
                                          e_array[0:i], e_l2_array[0:i], e_linf_array[0:i], e_min_array[0:i],
                                          maj_array[0:i], min_array[0:i])

if __name__ == '__main__':

    # ---------------------------------------------------------------------------------------------------------------------#
    # Define different parameters
    # ---------------------------------------------------------------------------------------------------------------------#

    # Define domain type
    unit_domain_tag = "unit-domain"
    l_shape_domain_tag = "l-shape-domain"
    pi_shape_domain_tag = "pi-shape-domain"
    circle_domain_tag = "circle-domain"

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

    # ---------------------------------------------------------------------------------------------------------------------#
    # Define a test form the dictionary of tests
    # ---------------------------------------------------------------------------------------------------------------------#
    tests = {1: tests.quadratic_polynomial_solution_1d_sigma_1, # checked
             2: tests.quadratic_polynomial_solution_1d_sigma_1000, # checked
             3: tests.quadratic_polynomial_solution_2d,
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
             44: tests.example_2_eps_1eminus3book_stynes_2d, # check if at least the indicator works fine
             45: tests.example_6_kleiss_tomar_2d, # doesn't coincide with work of Tomar and Kleiss
             46: tests.quadratic_polynomial_solution_2d_with_A_eps_eminus2,
             47: tests.example_1_sicom_paper,
             48: tests.neumann_bc_top_boundary_2d,
             49: tests.fokker_plank_example_1_2d, # to  test for fokker-plank paper
             50: tests.example_2_sicom_paper,
             51: tests.example_3_sicom_paper,
             52: tests.example_1_iga_majorant_test,
             53: tests.long_chen_example_1_2d, # eps = 1e-8
             54: tests.long_chen_example_2_2d, # eps = 1e-8
             55: tests.long_chen_example_3_2d, # eps = 1e-4
             56: tests.long_chen_example_1_eps_1minus1_2d, # eps = 1e-1 (to test different ranges to eps)
             57: tests.long_chen_example_1_eps_1minus2_2d, # eps = 1e-2
             58: tests.long_chen_example_1_eps_1minus3_2d,
             59: tests.long_chen_example_1_eps_1minus4_2d,
             60: tests.long_chen_example_1_eps_1minus5_2d,
             61: tests.long_chen_example_2_eps_1minus1_2d,
             62: tests.long_chen_example_2_eps_1minus2_2d,
             63: tests.long_chen_example_2_eps_1minus3_2d,
             64: tests.long_chen_example_2_eps_1minus4_2d,
             65: tests.long_chen_example_3_eps_eminus1_2d,
             66: tests.long_chen_example_3_eps_eminus2_2d,
             67: tests.long_chen_example_3_eps_eminus3_2d,
             68: tests.example_Neumann_BC_1d,
             69: tests.example_Neumann_BC_2d,
             70: tests.fokker_plank_example_2_2d
             }
    # ---------------------------------------------------------------------------------------------------------------------#
    # Set the number of the test and call for the problem data
    # ---------------------------------------------------------------------------------------------------------------------#
    test_num = 70
    
    u_expr, grad_u_expr, f_expr, \
    A_expr, lambda_1, lambda_expr, a_expr, uD_expr, uN_expr, \
    domain, dim, solution_tag, material_tag, pde_tag = tests[test_num]()

    # ---------------------------------------------------------------------------------------------------------------------#
    # Set the number of the test and call for the problem data
    # ---------------------------------------------------------------------------------------------------------------------#
    test_params = {'refinement_tag': adaptive_tag,
                   'nx0': 16, 'nx1': 16, 'nx2': 16, 'res':  16,
                   'marking_tag': bulk_tag,
                   'percentage_value': 0.4,
                   'refinement_criteria_tag': majorant_tag,
                   'number_of_refinements': 8,
                   'accuracy_level': 1e-3,
                   'v_approx_order': 1,
                   'flux_approx_order': 2,
                   'v_exact_approx_order': 4,
                   'solution_tag': solution_tag,
                   'material_tag':material_tag,
                   'pde_tag': pde_tag,
                   'majorant_optimization_iterations': 3,
                   'error_estimates': True,
                   'MAJORANT_OPTIMIZE': True,
                   'PLOT': True,
                   'STABILIZE': True,
                   'LOG_TO_FILE': False} # option that controls whether to log results into the file (True)
                                         # or to console (False)

    # ---------------------------------------------------------------------------------------------------------------------#
    # Run the test
    # ---------------------------------------------------------------------------------------------------------------------#
    test = TestEstimates(test_num,
                         u_expr, grad_u_expr, f_expr, A_expr, lambda_1, lambda_expr, a_expr, uD_expr, uN_expr,
                         domain, dim)
    test.test_estimates(test_params)




