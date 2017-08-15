__author__ = 'svetlana'

from dolfin.cpp.mesh import UnitIntervalMesh, cells, CellFunction, Point, UnitSquareMesh, UnitCubeMesh, Mesh
from dolfin import *
from dolfin.cpp.common import info, tic, toc, set_log_level
import time
from mshr import *


# Import functionality from the folder 'lib'
import estimates, errors, postprocess, integrators, problem
import tests


parameters["form_compiler"]["cpp_optimize"] = True
parameters["allow_extrapolation"] = True
#parameters["form_compiler_parameters"]["quadrature_degree"] = 3
set_log_level(PROGRESS)

def execute_refinement_strategy(problem_params, mesh, e_distr, m_distr, error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_max_k, h_min_k, t_T, k, nt,
                                project_path, results_folder, init_data_info):

    if problem_params['refinement_tag'] == "adaptive":
        # Define the distribution (error or majorant) upon which the refinement is based
        if problem_params['refinement_criteria_tag'] == 'error':
            distr = e_distr
        elif problem_params['refinement_criteria_tag'] == 'majorant':
            distr = m_distr

        if problem_params['marking_tag'] == 'average':
            marking = estimates.averaged_marking(mesh, distr)
            refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking'

        elif problem_params['marking_tag'] == 'predefined':
            marking = estimates.predefined_amount_of_elements_marking(mesh, distr, problem_params['percentage_value'])
            refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking' + '-proc-%d' % (problem_params['percentage_value'] * 100)

        elif problem_params['marking_tag'] == 'bulk':
            marking = estimates.bulk_marking(mesh, distr, problem_params['percentage_value'])
            refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking' + '-bulk-%d' % (problem_params['percentage_value'] * 100)


        mesh = refine(mesh, marking)

    elif problem_params['refinement_tag'] == "uniform":

        cell_markers = CellFunction("bool", mesh)
        cell_markers.set_all(True)
        mesh = refine(mesh, cell_markers)
        refinement_tag = problem_params['refinement_tag']

    num_cells = mesh.num_cells()
    num_vertices = mesh.num_vertices()

    postprocess.document_results_with_refinement(problem_params, mesh, e_distr, m_distr, error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_max_k, h_min_k,
                                                     t_T, k, nt, refinement_tag, project_path, results_folder, init_data_info, dim)

    return mesh
def solve_problem_nd_t(mesh, V, VV, V_exact, VV_exact, H_div,
                       f, sq_f,
                       phi, grad_phi,
                       u_e, grad_u_e, sq_grad_u_e,
                       u_0, u0_boundary,
                       dim, t_T, domain, C_FD,
                       nx0, nx1, nx2, nt,
                       delta, test_num,
                       project_path, results_folder, problem_params):

    # Define the value of the time-step
    tau = float(t_T / nt)

    # Initialize times
    t_k = 0         # t_k
    t_kp1 = tau     # t_{k + 1}
    t_km1 = -tau    # t_{k - 1}

    # Initialize the time integration scheme
    order = 4  # time integration order 2, 3, 4
    quadrature = integrators.TimeIntegrator(order, t_k, t_kp1)

    # Allocate the space for the arrays
    ed_k, et_k, e_incr_k, error_k, m0_k, md_k, mf_k, maj_incr_k, majorant_k, beta_k, i_eff_maj_k = \
        estimates.allocate_space_for_error_and_majorant(nt)
    h_max_k = postprocess.allocate_array(nt)
    h_min_k = postprocess.allocate_array(nt)
    maj_II_incr_k = postprocess.allocate_array(nt)
    error_d_k = postprocess.allocate_array(nt)
    majorant_d_k = postprocess.allocate_array(nt)

    vd_k = postprocess.allocate_array(nt)
    vt_k = postprocess.allocate_array(nt)
    v_norm_incr_k = postprocess.allocate_array(nt)
    v_norm_k = postprocess.allocate_array(nt)

    rel_error_k = postprocess.allocate_array(nt)
    rel_majorant_k = postprocess.allocate_array(nt)


    primal_problem_time = postprocess.allocate_array(nt)
    majorant_minimization_time = postprocess.allocate_array(nt)
    majorant_reconstruction_time = postprocess.allocate_array(nt)

    # Initialize the approximation and flux on the t = t_k
    v_k = interpolate(phi, V)
    y_k = interpolate(grad_phi, H_div)
    v_k1 = Function(V)

    # Initialize the counter
    k = 0

    while k+1 <= nt:

        t0_problem = tic()

        if k == 0 or problem_params['solving_strategy_tag'] == "with-refinement":
            # Define unknown function and test function
            u = TrialFunction(V)  # unknown function
            v = TestFunction(V)  # test function for the variational form

            # Define stiffness and mass matrices based on the functional space
            K, M = problem.construct_stiffness_and_mass_matrices(u, v, dim)
            # Define BC
            bc = DirichletBC(V, u_0, u0_boundary)

        # Update the time integration for the quadratures
        print "[%f, %f]:" % (t_k, t_kp1)
        quadrature.update(t_k, t_kp1)


        # Update right hand side f_{k}
        f.t = t_k
        f_k = interpolate(f, V)

        # Update right hand side f_{k + 1}
        f.t = t_kp1
        f_k1 = interpolate(f, V)

        # Update boundary condition
        u_0.t = t_kp1

        # Update the exact solution at t_{k + 1}
        u_e.t = t_kp1
        u_k1 = interpolate(u_e, V_exact)

        v_k1 = problem.get_next_step_solution(v_k, v_k1, K, tau, M, f_k, f_k1, bc, problem_params['time_discretization_tag'])

        #u_0.t = t_kp1
        #v_k1 = interpolate(u_e, V)

        if (problem_params['solving_strategy_tag'] == 'without-refinement' and k == 0) or \
            problem_params['solving_strategy_tag'] == 'with-refinement':
            postprocess.plot_solution(v_k1, mesh, dim, project_path + results_folder + 'u-%d' % (k+1))

        grad_v_k = problem.Grad(v_k, dim)
        grad_v_k1 = problem.Grad(v_k1, dim)

        # Define y_1 as a derivative with respect to x_i = x_0 = x[0]
        y_k1 = project(grad_v_k1, H_div)

        # Check the exact flux (test)
        #grad_phi.t = t_kp1
        #y_k1 = project(grad_phi, H_div)

        v_norm_incr_k, vd_k, vt_k = errors.v_norm_nd_t(k, v_norm_incr_k, vd_k, vt_k,
                                                       v_k1, grad_v_k, grad_v_k1,
                                                       V_exact, VV_exact,
                                                       mesh, dim, tau, delta)


        # Calculate error components
        e_incr_k, ed_k, et_k = errors.error_norm_nd_t(k, e_incr_k, ed_k, et_k,
                                                      grad_u_e, sq_grad_u_e, quadrature,
                                                      u_k1, v_k1, grad_v_k, grad_v_k1,
                                                      V_exact, VV_exact,
                                                      mesh, dim, tau, delta)
        # Calculate majorant components
        maj_incr_k, m0_k, md_k, mf_k, beta_k, y_k1 = \
            estimates.majorant_nd_t(k, maj_incr_k, m0_k, md_k, mf_k, beta_k, e_incr_k[k], ed_k[k], et_k[k],
                                    f, sq_f, phi, quadrature,
                                    v_k, v_k1, grad_v_k, grad_v_k1, V_exact,
                                    y_k, y_k1, H_div,
                                    mesh, tau, delta, dim, domain, C_FD,
                                    majorant_reconstruction_time, majorant_minimization_time)

        e_distr, m_distr, ed, md = estimates.error_majorant_distribution_nd(mesh, dim,
                                                                            grad_u_e, sq_grad_u_e, quadrature,
                                                                            grad_v_k, grad_v_k1, V_exact, VV_exact,
                                                                            y_k, y_k1, tau)

        # Document indicator perpormance
        print "h_max = ", mesh.hmax()
        h_max_k[k] = mesh.hmax()
        h_min_k[k] = mesh.hmin()

        init_data_info = postprocess.construct_result_tag(test_num, nx0, nx1, nx2, nt, k+1)

        if k%1 == 0 and problem_params['solving_strategy_tag'] == "with-refinement":
            #postprocess.plot_histogram(mesh, e_distr, m_distr,
            #                           project_path + results_folder + 'e-maj-distr-hist-' + problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking-' + results_info)
            postprocess.save_distr_to_mat_file(e_distr, m_distr,
                                               project_path + results_folder + 'e-maj-distr-hist-' + problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking-' + init_data_info)
        # Update the overall error from time interval [0, t_k]
        error_k, majorant_k, i_eff_maj_k, rel_error_k, rel_majorant_k = \
            estimates.add_increment_of_error_and_majorant(k, e_incr_k, maj_incr_k, error_k, et_k, majorant_k, i_eff_maj_k,
                                                          v_norm_incr_k, v_norm_k, vt_k, rel_error_k, rel_majorant_k)
        postprocess.output_time_layer_result_error_and_majorant(t_k, t_kp1, rel_error_k[k], rel_majorant_k[k], i_eff_maj_k[k])


        # Define the solving strategy (refine or not refine)
        if problem_params['solving_strategy_tag'] == "with-refinement":
            if problem_params['refinement_criteria_tag'] == "majorant":
                postprocess.plot_carpet_2d(mesh, ed, project_path + results_folder + 'carpet-error' + init_data_info)
                postprocess.plot_carpet_2d(mesh, md, project_path + results_folder + 'carpet-majorant' + init_data_info)


            # Refine mesh
            mesh = execute_refinement_strategy(problem_params, mesh, e_distr, m_distr,
                                               error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_max_k, h_min_k, t_T, k, nt,
                                               project_path, results_folder, init_data_info)
            num_cells = mesh.num_cells()
            num_vertices = mesh.num_vertices()
            #results_info = postprocess.construct_result_tag(test_num, nt, nx0, nx1, nx2, k+1, num_cells, num_vertices, dim)

            # Update functional spaces, BC, and stiffness/mass matrices
            V, VV, V_exact, VV_exact, H_div = problem.functional_spaces(mesh, problem_params["v_approx_order"],
                                                                        problem_params["flux_approx_order"], dim)
            # Update functions
            v_k.assign(interpolate(v_k1, V))
            y_k.assign(interpolate(y_k1, H_div))
            #y_k.assign(project(y_k1, H_div))

            v_k1 = Function(V)

        else:
            # Document results
            postprocess.document_results_without_refinement(mesh, e_distr, m_distr, error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_max_k, h_min_k, t_T, k, nt,
                                                            project_path, results_folder, init_data_info)
            # Update functions
            v_k.assign(v_k1)
            y_k.assign(y_k1)

        # Update time layer
        t_kp1 += tau
        t_km1 += tau
        t_k += tau
        k += 1

#------------------------------------------------------------
# nd x time problem
#------------------------------------------------------------

class TestMajorant():
    def __init__(self, test_num, u_expr, grad_u_expr, f_expr, u0_boundary, domain, dim, T):
        self.u_expr = u_expr
        self.grad_u_expr = grad_u_expr
        self.f_expr = f_expr
        self.u0_boundary = u0_boundary
        self.dim = dim
        self.domain = domain
        self.T = T

        self.test_num = test_num

    def convert_problem_data_to_expressions(self):

        sq_grad_u_expr = ""

        for i in range(0, self.dim):
            grad_u_x = self.grad_u_expr[i]
            index = grad_u_x.find(':')
            if index == -1:
                sq_grad_u_expr += 'pow( ' + grad_u_x[0:len(grad_u_x)] + ' , 2) + '
            else:
                sq_grad_u_expr += 'pow( ' + grad_u_x[index+1:len(grad_u_x)] + ' , 2) + '

        # Add the condition on the singularity before the squared
        if index != -1:
            condition = grad_u_x[0:index+1]
            sq_grad_u_expr = condition + '(' + sq_grad_u_expr[0 : len(sq_grad_u_expr) - 3] + ')'
        else:
            # Remove the plus sign in the end
            sq_grad_u_expr = sq_grad_u_expr[0:len(sq_grad_u_expr)-3]

        sq_f_expr = 'pow((' + self.f_expr + '), 2)'

        # Define exact solution expression
        u_e = Expression(self.u_expr, t=0)
        if self.dim == 1:
            grad_u_e = Expression(self.grad_u_expr[0], t=0)
        elif self.dim == 2:
            grad_u_e = Expression((self.grad_u_expr[0], self.grad_u_expr[1]), t=0)
        elif self.dim == 3:
            grad_u_e = Expression((self.grad_u_expr[0], self.grad_u_expr[1], self.grad_u_expr[2]), t=0)

        sq_grad_u_e = Expression(sq_grad_u_expr, t=0)

        # Define initial state
        phi = u_e
        grad_phi = grad_u_e

        # Define right-hand side
        f = Expression(self.f_expr, t=0)
        sq_f = Expression(sq_f_expr, t=0)

        return f, sq_f, phi, grad_phi, u_e, grad_u_e, sq_grad_u_e

    def test_majorant(self, problem_params):

        # Defnie the file path and project path
        project_path = postprocess.get_project_path()

        # Check if exist and create folder for the results
        results_folder_name = postprocess.construct_result_folder_name(self.dim, test_num, problem_params)

        postprocess.create_results_folder(results_folder_name)

        # Define mesh on the space domain
        if self.domain == "unit-domain":
            # 4 min
            nx0 = 4
            nx1 = 4
            nx2 = 4

            if self.dim == 1:
                # Create mesh on the unit interval
               mesh = UnitIntervalMesh(nx0)
            elif self.dim == 2:
                # Create mesh on the unit square
                mesh = UnitSquareMesh(nx0, nx1)
            elif self.dim == 3:
                # Create mesh on the unit cube
                mesh = UnitCubeMesh(nx0, nx1, nx2)

        elif self.domain == "l-shape-domain":
            nx0 = 0
            nx1 = 0
            nx2 = 0

            if dim == 2:

                data_path = '/data/%dd/test-%d/' % (self.dim, test_num)
                #file_name = 'l-shape-3d-mmg3d.xml'
                #file_name = 'l-shape-diam-2sqrt2-33-vertices.xml'
                file_name = 'l-shape-3.xml'
                #file_name = 'l-shape-15-verts.xml'
                #file_name = 'l-shape-adapt.xml'

                #file_name = 'l-shape-3d-twice-refined.xml'
                #file_name = 'l-shape-3d-once-refined.xml'
                #plot(mesh, interactive=True)

                # Load mesh from the xml-file
                mesh = Mesh(project_path + data_path + file_name)
                plot(mesh, interactive=True, mode="glyphs", scale=2.0)


                '''
                res = 32
                big_square = Rectangle(Point(-1.0, -1.0), Point(1.0, 1.0))
                small_square = Rectangle(Point(0.0, -1.0), Point(1.0, 0.0))
                mesh = generate_mesh(big_square - small_square, res)
                plot(mesh, interactive=True, mode="glyphs", scale=2.0)
                '''

                # Plot initial mesh
            elif dim == 3:
                mesh = generate_mesh(Box(Point(0, 0, 0), Point(1, 1, 1)) +
                                     Box(Point(0, 0, 0), Point(-1, 1, 1)) +
                                     Box(Point(0, 0, 0), Point(-1, -1, 1)), 4)
                plot(mesh, interactive=True)

        elif self.domain == "pi-shape-domain":
            nx0 = 0
            nx1 = 0
            nx2 = 0

            data_path = '/data/%dd/test-%d/' % (self.dim, test_num)
            file_name = 'pi-shape-2.xml'

            # Load mesh from the xml-file
            mesh = Mesh(project_path + data_path + file_name)
            plot(mesh, interactive=True)

        elif self.domain == "circle-domain":

            mesh = generate_mesh(Circle(dolfin.Point(0, 0), 1), 32)
            plot(mesh, interactive=True)


        if dim == 2:
            postprocess.plot_mesh(mesh, project_path + results_folder_name + 'initial-mesh')
        elif dim == 3:
            postprocess.plot_mesh_3d(mesh, project_path + results_folder_name + 'initial-mesh')

        # Mesh parameters in time dimention
        nt = 10 # second iteration
        #nt = 80 # second-third iteration
        #nt = 160 # 3rd
        #nt = 320 #
        #nt = 1280 #

        # Define mesh and functional spaces based on the mesh
        V, VV, V_exact, VV_exact, H_div = problem.functional_spaces(mesh,
                                                                    problem_params["v_approx_order"],
                                                                    problem_params["flux_approx_order"],
                                                                    self.dim)

        # Define problem data functions
        f, sq_f, phi, grad_phi, u_e, grad_u_e, sq_grad_u_e = self.convert_problem_data_to_expressions()

        # Dirichlet boundary condition (x = 0 or x = 1)
        u_0 = phi

        # Majorant parameter
        delta = 1.0

        # Output to console the problem data
        postprocess.output_problem_characteristics(self.test_num, self.u_expr, self.grad_u_expr, self.f_expr,
                                                   self.T, self.domain, self.dim,
                                                   nt, mesh.num_cells(), mesh.num_vertices(),
                                                   problem_params["v_approx_order"],
                                                   problem_params["flux_approx_order"],
                                                   delta)

        # Calculate the estimate for Freidrichs constant based on the domain
        C_FD = problem.calculate_CF_of_domain(domain, self.dim)

        t0_problem = time.clock()

        solve_problem_nd_t(mesh, V, VV, V_exact, VV_exact, H_div,
                           f, sq_f,
                           phi, grad_phi,
                           u_e, grad_u_e, sq_grad_u_e,
                           u_0, self.u0_boundary,
                           self.dim, self.T, self.domain, C_FD,
                           nx0, nx1, nx2, nt,
                           delta,
                           self.test_num,
                           project_path, results_folder_name, 
                           problem_params)

        t1_problem = time.clock()
        print "time = ", t1_problem - t0_problem

# Dictionary of tests
tests = {1: tests.polynomial_solution_3d_t, # unit domain, u = (1 - x[0])*x[0]*(1 - x[1])*x[1]*(1 - x[2])*x[2]*(t*t + t + 1)
         2: tests.solution_with_singularities_3d_t, # l-shape domain in 3d, u with singularities
         3: tests.solution_without_singularities_3d_t, # l-shape domain in 3d, u without singularities
         4: tests.solution_with_singularities_2d_t,
         5: tests.polynomial_solution_2d_t,
         6: tests.solution_without_singularities_2d_t, # pi-shape smail in 2d
         7: tests.polynomial_in_space_trigonometric_in_time_solution_1d_t,
         8: tests.trigonometric_drastic_changing_in_time_2d_t,
         9: tests.polynomial_in_space_exponential_in_time_1d_t,
         10: tests.polynomial_in_space_linear_in_time_1d_t,
         11: tests.solution_with_singularities_sint_t_2d_t}
# Set the number of the test and call for the problem data
test_num = 4
u_expr, grad_u_expr, f_expr, domain, dim, T = tests[test_num]()

# Define domain type
unit_domain_tag = "unit-domain"
l_shape_domain_tag = "l-shape-domain"
pi_shape_domain_tag = "pi-shape-domain"

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

# Define solving strategy
with_refinement_tag = "with-refinement"
without_refinement_tag = "without-refinement"

# Define solving strategy
material_changing_tag = "material-changing"
material_constant_tag = "material-constant"

# Define solving strategy
explicit_tag = "explicit"
implicit_tag = "implicit"

# Functional spaces parameters
v_deg = 1
y_deg = 2

# Init problem parameters
problem_params = dict(refinement_tag=adaptive_tag,
                      marking_tag=bulk_tag,
                      percentage_value=0.3,
                      refinement_criteria_tag=majorant_tag,
                      solving_strategy_tag=with_refinement_tag,
                      v_approx_order=v_deg,
                      flux_approx_order=y_deg,
                      time_discretization_tag=implicit_tag)

# Define Dirichlet boundary condition
def u0_boundary(x, on_boundary):
    return on_boundary

test = TestMajorant(test_num, u_expr, grad_u_expr, f_expr, u0_boundary, domain, dim, T)
test.test_majorant(problem_params)




