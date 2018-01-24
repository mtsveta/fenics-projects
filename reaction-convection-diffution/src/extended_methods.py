__author__ = 'svetlana'

from dolfin import *
from dolfin.cpp.mesh import cells
import postprocess, integration
import numpy


def calculate_majorant(u, y, beta, f, lmbd, mesh, dim, C_FD):

    # Define mu parameter
    mu = 1.0

    # Define residual of majorant
    r_d = (Grad(u, dim) - y)
    r_f = (Div(y, dim) + f - lmbd * u)

    # Define majorant components
    m_d = assemble(inner(r_d, r_d) * dx(domain=mesh))
    m_f_mu_opt = assemble((mu ** 2) / lmbd * inner(r_f, r_f) * dx(domain=mesh))
    m_f_one_minus_mu_opt = assemble(((1 - mu) ** 2) * inner(r_f, r_f) * dx(domain=mesh))
    m_f = m_f_mu_opt + (1.0 + 1.0 / beta) * (C_FD ** 2) * m_f_one_minus_mu_opt

    # Calculate majorant based on the parameter value
    if m_f <= DOLFIN_EPS:
        maj = m_d
    else:
       if m_d <= DOLFIN_EPS:
           maj = m_f
       else:
           # Calculate the optimal value for beta parameter
           beta = C_FD * sqrt(m_f / m_d)
           maj = (1.0 + beta) * m_d + m_f

    return maj, m_d, m_f_one_minus_mu_opt, beta

def error_majorant_distribution_nd(mesh, domain, dim,
                                   u_e, grad_u_e, u, V_exact, y, f, lmbd, beta):

    cell_num = mesh.num_cells()
    #e_distr = postprocess.allocate_array(cell_num)
    #m_distr = postprocess.allocate_array(cell_num)
    ed_distr = postprocess.allocate_array(cell_num)
    md_distr = postprocess.allocate_array(cell_num)

    # Define error and residuals
    r_d = abs(Grad(u, dim) - y)
    r_f = abs(f + Div(y, dim) - lmbd * u)

    #C_FD = calculate_CF_of_domain(domain, dim)

    u_ve = interpolate(u, V_exact)
    u_exact_ve = interpolate(u_e, V_exact)
    e = u_ve - u_exact_ve

    #w_opt = (C_FD ** 2) * (1 + beta) / (beta + (C_FD ** 2) * (1 + beta) * lmbd)

    # Project UFL forms on the high-order functional space to obtain functions
    #ed = project(sqrt(inner(Grad(e, dim), Grad(e, dim))), V_exact)
    #md = project(sqrt(inner(r_d, r_d)), V_exact)
    ed = project(inner(Grad(e, dim), Grad(e, dim)), V_exact)
    md = project(inner(r_d, r_d), V_exact)

    #ed_test = assemble(sqrt(inner(Grad(e, dim), Grad(e, dim))) * dx)
    #md_test = assemble(sqrt(inner(r_d, r_d)) * dx)
    #e = project(sqrt(inner(Grad(e, dim), Grad(e, dim)) + lmbd * inner(e, e)), V_exact)
    #m = project(sqrt((1 + beta) * inner(r_d, r_d) + w_opt * inner(r_f, r_f)), V_exact)

    dofmap = V_exact.dofmap()

    scheme_order = 4
    gauss = integration.SpaceIntegrator(scheme_order, dim)

    for c in cells(mesh):
        verts = dofmap.tabulate_coordinates(c)
        x_n = verts[0:dim+1, :]
        x_n_transpose = x_n.transpose()

        matrix = postprocess.allocate_array_2d(dim+1, dim+1)
        matrix[0:dim, 0:dim+1] = x_n_transpose
        matrix[dim, 0:dim+1] = numpy.array([1.0 for i in range(dim + 1)])
        meas = abs(numpy.linalg.det(matrix))

        ed_distr[c.index()] = gauss.integrate(ed, meas, x_n_transpose)
        md_distr[c.index()] = gauss.integrate(md, meas, x_n_transpose)

        #e_distr[c.index()] = gauss.integrate(e, meas, x_n_transpose)
        #m_distr[c.index()] = gauss.integrate(m, meas, x_n_transpose)

    print 'sum of m_cells', numpy.sum(md_distr)
    print 'sum of e_cells', numpy.sum(ed_distr)

    return ed_distr, md_distr


def majorant_distribution_DG0(mesh, V, V_exact, f, u0, u, u_exact, y, beta, C_FD, dim):

    # Define the functional space used for distribution over cells
    DG0 = FunctionSpace(mesh, "DG", 0)
    w = TestFunction(DG0)
    m_d_DG0 = Function(DG0)
    m_df_DG0 = Function(DG0)

    r_d = abs(Grad(u, dim) - y)
    r_f = abs(f + Div(y, dim))
    # Define variational forms of dual component of majorant and the whole functional
    m_d_var = w * sqrt(inner(r_d, r_d)) * dx
    m_df_var = w * sqrt((1 + beta) * inner(r_d, r_d) + C_FD**2 * (1 + 1 / beta) * inner(r_f, r_f)) * dx

    # Assemble the variation form and dumping obtained vector into function from DG0
    assemble(m_d_var, tensor=m_d_DG0.vector())
    assemble(m_df_var, tensor=m_df_DG0.vector())

    # Get the distribution of the indicator and majorant over cells
    m_d_cells = m_d_DG0.vector().array()
    m_df_cells = m_df_DG0.vector().array()

    #print "m_d_DG0 total", numpy.sum(m_d_cells)
    #print "m_df_DG0 total", numpy.sum(m_df_cells)

    return m_d_cells, m_df_cells



# goal functional - example 1 (green book)
def J(w, grad_e, norm_grad_e, dim):
    return inner(Grad(w, dim), grad_e) / norm_grad_e

def solve_dual_problem(V_star, f, u0, u0_boundary, u, u_e, dim):

    # Define dual variational problem
    z = TrialFunction(V_star)
    w = TestFunction(V_star)

    # Define the system
    a_star = inner(Grad(z, dim), Grad(w, dim)) * dx(domain=mesh)
    grad_e = Grad(u_e, dim) - Grad(u, dim)
    norm_grad_e = sqrt(assemble(inner(grad_e, grad_e) * dx(domain = mesh)) )
    L_star = J(w, grad_e, norm_grad_e, dim)

    z = Function(V_star)
    bc_star = DirichletBC(V_star, u0, u0_boundary)
    solve(a_star == L_star, z, bc_star)

    z_exact = project( (u - u_e) / norm_grad_e, V_star)
    z_error = assemble((z - z_exact)**2 * dx(domain=mesh))

    print "z error", z_error

    return z

def compare_error_indicators(mesh, V, V_star, f, u0, u0_boundary, u, u_e, grad_u_e, y, beta, test_num, tag, dim):

    z = solve_dual_problem(V_star, f, u0, u0_boundary, u, u_e)
    Pz = project(z, V)

    norm_grad_e = sqrt(assemble(inner(grad_u_e - Grad(u, dim), grad_u_e - Grad(u, dim)) * dx(domain = mesh)))
    h = mesh.hmax()
    n = FacetNormal(mesh)
    DG0 = FunctionSpace(mesh, "DG", 0)
    w = TestFunction(DG0)

    r_h = jump(Grad(u, dim), n)
    R_h = div(Grad(u, dim)) + f

    cell_num = mesh.num_cells()

    e = abs(u_e - u)
    r_d = abs(Grad(u, dim) - y)
    r_f = abs(f + Div(y, dim))
    eta_var = w * (h * R_h)**2 * dx + avg(w) * avg(h) * r_h**2 * dS
    E_DWR_var = w * inner(R_h, z - Pz) * dx - 0.5 * avg(w) * inner(r_h, avg(z - Pz)) * dS
    J_e_var = w * J(w, e, norm_grad_e) * dx
    m_d_var = w * sqrt(inner(r_d, r_d)) * dx
    #C_FD = 1 / (sqrt(3) * DOLFIN_PI)
    height = 1
    width = 2
    length = 2
    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height**2 + 1.0 / width**2 + 1.0 / length**2)
    m_df_var = w * sqrt((1 + beta) * inner(r_d, r_d) + C_FD**2 * (1 + 1 / beta) * inner(r_f, r_f)) * dx

    #M_e = w * (u - u_exact) * dx

    eta_DG0 = Function(DG0)
    E_DWR_DG0 = Function(DG0)
    J_e_DG0 = Function(DG0)
    m_d_DG0 = Function(DG0)
    m_df_DG0 = Function(DG0)

    assemble(eta_var, tensor=eta_DG0.vector())
    assemble(E_DWR_var, tensor=E_DWR_DG0.vector())
    assemble(J_e_var, tensor=J_e_DG0.vector())
    assemble(m_d_var, tensor=m_d_DG0.vector())
    assemble(m_df_var, tensor=m_df_DG0.vector())

    eta_distr = eta_DG0.vector().array()
    E_DWR_distr = E_DWR_DG0.vector().array()
    J_e_distr = J_e_DG0.vector().array()
    m_d_distr = m_d_DG0.vector().array()
    m_df_distr = m_df_DG0.vector().array()

    eta_DG0_total = numpy.sum(eta_distr)
    E_DWR_DG0_total = numpy.sum(E_DWR_distr)
    J_e_DG0_total = numpy.sum(J_e_distr)
    m_d_DG0_total = numpy.sum(m_d_distr)
    m_df_DG0_total = numpy.sum(m_df_distr)

    print "eta_DG0 total:", eta_DG0_total
    print "E_DWR_DG0 total", E_DWR_DG0_total
    print "J_e_DG0 total", J_e_DG0_total
    print "m_d_DG0 total", m_d_DG0_total
    print "m_df_DG0 total", m_df_DG0_total

    #plot_histogram(1, 0, cell_num, 0, V, mesh, J_e_DG0.vector().array(), E_DWR_DG0.vector().array(), 'dwr-')
    #plot_histogram(1, 0, cell_num, 0, V, mesh, J_e_DG0.vector().array(), eta_DG0.vector().array(), 'residual-')
    #plot_histogram(1, 0, cell_num, 0, V, mesh, J_e_DG0.vector().array(), m_d_DG0.vector().array(), 'maj-')
    '''
    save_distr_to_mat_file(test_num, cell_num, V,
                           eta_distr,
                           E_DWR_distr,
                           J_e_distr,
                           m_d_distr,
                           m_df_distr,
                           tag + '-dg0-')

    '''
    #return eta_distr, E_DWR_distr, J_e_distr, m_d_distr, m_df_distr


    '''
    cell_domains = CellFunction('size_t', mesh)
    facet_domains = FacetFunction('size_t', mesh, 0)

    eta_omega_cell = CellFunction('double', mesh)
    E_DWR_cell = CellFunction('double', mesh)
    J_e_cell = CellFunction('double', mesh)
    M_cell = CellFunction('double', mesh)

    eta_omega = 0
    E_DWR = 0
    J_e = 0
    M = 0

    cell_domains.set_all(0)
    facet_domains.set_all(0)

    for cell in cells(mesh):

        indx = cell.index() + 1
        cell_domains[cell] = indx
        h = cell.diameter()

        # Assemble components of error indicator over cell
        R_K = assemble((h*R_h)**2 * dx(indx), cell_domains=cell_domains)
        Z_K = assemble((1/h*(z - Pz))**2 * dx(indx), cell_domains=cell_domains)
        RZ_K = assemble(inner(R_h, z - Pz) * dx(indx), cell_domains=cell_domains)
        J_e_K = assemble(inner(grad(e), grad(e)) / inv_const * dx(indx), cell_domains=cell_domains)
        m_d_K = assemble(sqrt(inner(grad(u) - y, grad(u) - y)) * dx(indx), cell_domains=cell_domains)

        for facet in facets(cell):
            if facet_domains[facet] == 0 and not(facet.exterior()):
                facet_domains[facet] = indx

        # Assemble components of error indicator over facets of the cell
        r_K = assemble(avg(h) * r_h**2 * dS(indx), interior_facet_domains=facet_domains)
        z_K = assemble(1/ avg(h) * avg(z - Pz)**2 * dS(indx), interior_facet_domains=facet_domains)
        rz_K = assemble(inner(r_h, avg(z - Pz)) * dS(indx), interior_facet_domains=facet_domains)

        rho_K = sqrt(R_K + r_K)
        omega_K = sqrt(Z_K + z_K)
        e_dwr_K = RZ_K - 0.5 * rz_K
        #m_d_K = sqrt(m2_d_K)

        eta_omega_cell[cell] = rho_K * omega_K
        E_DWR_cell[cell] = e_dwr_K
        J_e_cell[cell] = J_e_K
        M_cell[cell] = m_d_K

        eta_omega += eta_omega_cell[cell]
        E_DWR += E_DWR_cell[cell]
        J_e += J_e_cell[cell]
        M += M_cell[cell]

    print "eta_omega over cells", eta_omega_cell.array()
    print "E_DWR over cells", E_DWR_cell.array()
    print "J_e over cells", J_e_cell.array()
    print "M over cells", M_cell.array()

    plot_histogram(1, 0, cell_num, 0, V, mesh, J_e_cell.array(), E_DWR_cell.array(), 'dwr-over-cells-')
    plot_histogram(1, 0, cell_num, 0, V, mesh, J_e_cell.array(), eta_omega_cell.array(), 'residual-over-cells-')
    plot_histogram(1, 0, cell_num, 0, V, mesh, J_e_cell.array(), M_cell.array(), 'maj-over-cells-')

    print "eta", eta_omega
    print "dwr", E_DWR
    print "er", J_e
    print "maj", M
    '''
    '''
    eta_var = w * (h * R_h)**2 * dx + avg(w) * avg(h) * r_h**2 * dS
    E_DWR_var = w * inner(R_h, z - Pz) * dx - 0.5 * avg(w) * inner(r_h, avg(z - Pz)) * dS
    J_e_var = w * inner(grad(e), grad(e)) / inv_const * dx

    eta_DG0 = Function(DG0)
    E_DWR_DG0 = Function(DG0)
    J_e_DG0 = Function(DG0)
    m_d_DG0 = Function(DG0)
    m_df_DG0 = Function(DG0)

    assemble(eta_var, tensor=eta_DG0.vector())
    assemble(E_DWR_var, tensor=E_DWR_DG0.vector())
    assemble(J_e_var, tensor=J_e_DG0.vector())
    assemble(m_d_var, tensor=m_d_DG0.vector())
    assemble(m_df_var, tensor=m_df_DG0.vector())

    eta_DG0_total = numpy.sum(eta_DG0.vector().array())
    E_DWR_DG0_total = numpy.sum(E_DWR_DG0.vector().array())
    J_e_DG0_total = numpy.sum(J_e_DG0.vector().array())
    m_d_DG0_total = numpy.sum(m_d_DG0.vector().array())
    m_df_DG0_total = numpy.sum(m_df_DG0.vector().array())

    print "eta_DG0 total:", eta_DG0_total
    print "E_DWR_DG0 total", E_DWR_DG0_total
    print "J_e_DG0 total", J_e_DG0_total
    print "m_d_DG0 total", m_d_DG0_total
    print "m_df_DG0 total", m_df_DG0_total

    #plot_histogram(1, 0, cell_num, 0, V, mesh, J_e_DG0.vector().array(), E_DWR_DG0.vector().array(), 'dwr-')
    #plot_histogram(1, 0, cell_num, 0, V, mesh, J_e_DG0.vector().array(), eta_DG0.vector().array(), 'residual-')
    #plot_histogram(1, 0, cell_num, 0, V, mesh, J_e_DG0.vector().array(), m_d_DG0.vector().array(), 'maj-')

    eta_distr = eta_DG0.vector().array()
    E_DWR_distr = E_DWR_DG0.vector().array()
    J_e_distr = J_e_DG0.vector().array()
    m_d_distr = m_d_DG0.vector().array()
    m_df_distr = m_df_DG0.vector().array()

    return eta_distr, E_DWR_distr, J_e_distr, m_d_distr, m_df_distr

    '''



#-------------------------------------------------------------------------------#
def solve_problem_nd_t(mesh, V, VV, V_exact, VV_exact, H_div,
                       f, lmbd, phi, grad_phi,
                       u_e, grad_u_e,
                       u_0, u0_boundary,
                       dim, domain, C_FD, nx0, nx1, nx2,
                       test_num,
                       results_folder,
                       problem_params):

    print " "
    print "%---------------------"
    print "% Solving the problem "
    print "%---------------------"
    print " "

    # Initialize the project path
    project_path = postprocess.get_project_path()

    if problem_params['solving_strategy_tag'] == "with-refinement":
        ref_num = 5
    else:
        ref_num = 1

    ref_num_array = postprocess.allocate_array(ref_num)
    h_max_array = postprocess.allocate_array(ref_num)

    e_array = postprocess.allocate_array(ref_num)

    maj_array = postprocess.allocate_array(ref_num)
    min_array = postprocess.allocate_array(ref_num)

    i_eff_maj_array = postprocess.allocate_array(ref_num)
    i_eff_min_array = postprocess.allocate_array(ref_num)

    for i in range(0, ref_num):

        # Compute approximate solution, its error and majorant
        u = problem.solve_poisson(V, f, lmbd, u_0, u0_boundary, mesh)

        # Calculate error
        error = estimates.error_norm(u, u_e, lmbd, V_exact, mesh, dim)

        # Calculate majorant
        y = project(estimates.Grad(u, dim), H_div)
        MAJORANT_OPTIMIZE = 1

        #y = project(grad_phi, H_div)
        #MAJORANT_OPTIMIZE = 0

        # Calculate majorant
        maj, y, beta, md, mf, majorant_reconstruction_time, majorant_minimization_time = \
            estimates.majorant_nd(u, y, H_div, f, lmbd, error, mesh, dim, domain, C_FD, MAJORANT_OPTIMIZE)
        i_eff_maj = sqrt(maj / error)

        # Construct error and majorant distribution
        e_distr, m_distr = estimates.error_majorant_distribution_nd(mesh, dim, u_e, u, V_exact, y)


        print "% -------------"
        print "% Save results "
        print "% -------------"

        num_cells = mesh.num_cells()
        num_vertices = mesh.num_vertices()

        # Construct the tag with problem information
        results_info = postprocess.construct_result_tag(test_num, i, nx0, nx1, nx2, num_cells, num_vertices)

        # Plot histogram with the error and majorant distribution
        postprocess.plot_histogram(mesh, e_distr, m_distr,
                                   project_path + results_folder + 'e-maj-distr-hist' + results_info)
        # Output the results of upper bound reconstructions
        estimates.output_result_error_and_majorant(error, maj, i_eff_maj)

        min, phi = estimates.minorant(u, mesh, V_exact, u_0, u0_boundary, f, lmbd, dim, error)
        i_eff_min = sqrt(min / error)
        estimates.output_result_minorant(min, min/error)

        # Update the arrays of data with respect to the refinement cycle
        ref_num_array[i] = i + 1
        e_array[i] = error
        maj_array[i] = maj
        min_array[i] = min
        i_eff_maj_array[i] = i_eff_maj
        i_eff_min_array[i] = i_eff_min
        h_max_array[i] = mesh.hmax()

        # Define the refinement strategy
        if problem_params['refinement_tag'] == "adaptive":
            # Define the distribution (error or majorant) upon which the refinement is based
            if problem_params['refinement_criteria_tag'] == 'error':
                distr = e_distr
            elif problem_params['refinement_criteria_tag'] == 'majorant':
                distr = m_distr
            #elif test_params['refinement_criteria_tag'] == 'e-dwr':
            #    distr = e_dwr_distr
            #elif test_params['refinement_criteria_tag'] == 'residual':
            #    distr = residual_distr

            # Run the marking procedure
            theta = 0.4    # percentage of the marked cells

            if problem_params['marking_tag'] == 'average':
                marking = estimates.averaged_marking(mesh, distr)
            elif problem_params['marking_tag'] == 'predefined':
                marking = estimates.predefined_amount_of_elements_marking(mesh, distr, theta)
            elif problem_params['marking_tag'] == 'bulk':
                marking = estimates.bulk_marking(mesh, distr, theta)

            # Refine mesh based on the marking criteria
            mesh = refine(mesh, marking)

            refinement_tag = problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking'
            refinement_tag = refinement_tag + '-theta-%d' % (theta * 100)

            # Plot result mesh
            postprocess.plot_mesh(mesh, project_path + results_folder + 'mesh-based-' + refinement_tag + results_info)

            # Save the results in mat-file
            postprocess.save_results_to_mat_file(e_distr, m_distr, error, maj, i_eff_maj, min, i_eff_min,
                                                 project_path + results_folder + 'results-on-mesh-based-' +
                                                 refinement_tag + results_info)
        elif problem_params['refinement_tag'] == "uniform":

            # Plot result mesh
            postprocess.plot_mesh(mesh, project_path + results_folder + 'mesh-uniform' + results_info)

            # Save the results in mat-file
            postprocess.save_results_to_mat_file(e_distr, m_distr, error, maj, i_eff_maj, min, i_eff_min,
                                                 project_path + results_folder + 'results-on-mesh-uniform' + results_info)

            # Refinement for the comparison with theta < 1.0
            #cell_markers = CellFunction("bool", mesh)
            #cell_markers.set_all(True)
            #mesh = refine(mesh, cell_markers)

            # Refinement for the optimal convergence
            mesh = refine(mesh)

        if i == ref_num - 1:

            # Define the refinement_tag
            if problem_params['refinement_tag'] == "adaptive":
                refinement_tag = problem_params['refinement_criteria_tag'] + '-' + \
                                 problem_params['marking_tag'] + '-marking'
            elif problem_params['refinement_tag'] == "uniform":
                refinement_tag = problem_params['refinement_tag']

            postprocess.plot_decay_error_majorant(problem_params["v_approx_order"], h_max_array, e_array, maj_array,
                                                  project_path + results_folder +
                                                  "error-maj-v-P%d-" % problem_params["v_approx_order"] +
                                                  refinement_tag + results_info)
            postprocess.save_decay_to_mat_file(h_max_array, e_array, maj_array,
                                               project_path + results_folder +
                                               "maj-decay-v-P%d" % problem_params["v_approx_order"])
            postprocess.plot_decay_majorant_v_of_deg(h_max_array, maj_array, problem_params["v_approx_order"],
                                                     project_path + results_folder +
                                                     "maj-decay-v-P%d" % problem_params["v_approx_order"] +
                                                     results_info)

        # Update functional spaces, BC, and stiffness/mass matrices
        V, VV, V_exact, VV_exact, H_div = problem.update_functional_space(mesh,
                                                                          problem_params["v_approx_order"],
                                                                          problem_params["flux_approx_order"],
                                                                          dim)

        #tag = 'E-DWR-'
        #tag = 'eta-'
        #tag = 'error-'
        #tag = 'm-d-'
        tag = 'm-df-'

        '''
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