__author__ = 'svetlana'

import scipy.io
import os, sys
import numpy

from dolfin.cpp.mesh import cells
from dolfin import *
import plotting

import matplotlib.pylab as plt
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib
from matplotlib import cm

# Memory allocation functions
def allocate_array(n): return numpy.array([0.0 for i in range(n)])
def allocate_array_2d(m, n): return numpy.array([[0.0 for j in range(n)] for i in range(m)])

# -------------------------------------------------------------------------------------------------------------------- #
# Data translation functions
# -------------------------------------------------------------------------------------------------------------------- #

# Constructing the name of the result folder based on the problem parameters
def construct_result_folder_name(dim, test_num, problem_params):
    results_folder_name = '/out/%dd/test-%d/v-P%d-y-RT%d/' \
                          % (dim, test_num,
                             problem_params["v_approx_order"],
                             problem_params["flux_approx_order"])
    if problem_params["refinement_tag"] == 'uniform':
        results_folder_name += problem_params["refinement_tag"]
    else:
        results_folder_name += problem_params["marking_tag"] + '-' + str(100 * problem_params["percentage_value"])

    if problem_params["STABILIZE"]:
        results_folder_name += '-stabilized/'
    else:
        results_folder_name += '/'

    return results_folder_name

# Function creating the tag for the result folder
def construct_result_tag(test_num, ref_num, nx0, nx1, nx2, num_cells, num_vertices):
    file_info = '-test-num-%d-ref-inter-%d-%d-%d-%d-num-cells-%d-num-vertices-%d' \
                % (test_num, ref_num, nx0, nx1, nx2, num_cells, num_vertices)
    return file_info

# Function returning the path of the project
def get_project_path():
    (project_path, src_folder) = os.path.split(os.path.abspath(os.path.dirname(__file__)))
    #(project_path, src_folder) = os.path.split(os.path.abspath(src_path))
    return project_path

# Function creating the result folder
def create_results_folder(directory, test_params):
    # Creat the full directory of the result folder
    full_directory = get_project_path() + directory

    # If derictory doesn't exist
    if not os.path.exists(full_directory):
        os.makedirs(full_directory)
    # Write into the loggin path
    if test_params["LOG_TO_FILE"] == True:
        sys.stdout = open(full_directory + "results.txt", "w+")

def create_decay_tag(test_params, project_path, results_folder, results_info):
    if test_params['refinement_tag'] == "adaptive":
        refinement_tag = test_params['refinement_criteria_tag'] + '-' + \
                         test_params['marking_tag'] + '-marking' + '-' + str(100*test_params["percentage_value"])
    elif test_params['refinement_tag'] == "uniform":
        refinement_tag = test_params['refinement_tag']

    decay_result_folder = project_path + results_folder + \
                          "error-maj-v-P%d-" % test_params["v_approx_order"] \
                          + refinement_tag + results_info
    return decay_result_folder

def document_errors_decay(dim, test_params, decay_result_folder,
                          dofs_array, h_min_array, e_array, e_l2_array, e_linf_array, e_min_array, maj_array, min_array):

    output_decay_error_summary(test_params, dim, p, dofs_array, h_min_array, e_array, e_l2_array, e_linf_array)

    # Output errors with error estimates
    if test_params["error_estimates"] == True:
        if test_params["PLOT"] == True:
            plot_decay_error_estimates_dofs(dim, float(test_params["v_approx_order"]),
                                            dofs_array, h_min_array, e_array, maj_array, min_array,
                                            decay_result_folder)
        output_decay_error_estimates_summary(test_params, dim, float(test_params["v_approx_order"]),
                                            dofs_array, h_min_array, e_array, maj_array, min_array)
        output_error_estimates_summary(dofs_array, e_array, maj_array, e_min_array, min_array)
    # Output errors
    if test_params["PLOT"] == True:
        plot_decay_error_dofs(dim, float(test_params["v_approx_order"]), dofs_array, h_min_array, e_array, e_l2_array, e_linf_array,
                              decay_result_folder)

    save_to_mat_file(dofs_array, 'dofs_array', decay_result_folder)
    save_to_mat_file(maj_array, 'maj_array', decay_result_folder)
    save_to_mat_file(e_array, 'e_array', decay_result_folder)


# Saving data into the mat-file
def save_to_mat_file(data, data_tag, result_path):
    scipy.io.savemat(result_path, {data_tag: data})


def save_results_to_mat_file(error_distr, maj_distr, error, majorant, i_eff_maj, result_path):
    scipy.io.savemat(result_path,
                     {'error_distr': error_distr, 'maj_distr': maj_distr,
                      'error': error, 'majorant': majorant, 'i_eff_maj': i_eff_maj})

# Saveing result mesh into xml file
def save_mesh_to_xml_file(mesh, file_path, file_tag):
    dolfin.File(file_path + 'mesh' + file_tag + '.xml') << mesh
# -------------------------------------------------------------------------------------------------------------------- #
# Plotting functions
# -------------------------------------------------------------------------------------------------------------------- #

# Function plotting distribution of the error an  majorant
def plot_histogram(mesh, e_distr, m_distr, result_path):

    num_cells = mesh.num_cells()
    cells = numpy.linspace(0, num_cells, num_cells)

    plt.figure()

    plt.bar(cells, e_distr, align="center", facecolor="red")
    plt.plot(cells, m_distr, 'g.-')
    plt.xlabel('cells')
    plt.ylabel('e distr, m distr')
    plt.xlim([cells[0], cells[num_cells - 1]])

    fig = plt.gcf()
    fig.savefig(result_path + '.eps', format="eps")
    plt.close('all')

# Function plotting distribution of the error an  majorant
def plot_histogram_smoothness(mesh, smoothness_distr, color, result_path):

    num_cells = mesh.num_cells()
    cells = numpy.linspace(1, num_cells, num_cells)

    plt.figure()
    plt.bar(cells, smoothness_distr, align="center", facecolor=color)
    plt.xlabel('cells')
    plt.ylabel('')
    plt.xlim([cells[0], cells[num_cells-1]])

    fig = plt.gcf()
    fig.savefig(result_path + '.eps', format="eps")
    plt.close('all')

# Plotting the mesh
def plot_mesh(mesh, result_path):
    plt.figure()
    plt.hold(True)

    # change the font
    matplotlib.rcParams.update({'font.size': 10,
                                'font.family': 'serif'})

    n_vert = mesh.num_vertices()
    n_cells = mesh.num_cells()
    d = mesh.geometry().dim()

    plotting.plot(mesh)
    # Create the triangulation
    # mesh_coordinates = mesh.coordinates().reshape((n_vert, d))
    # triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    # triangulation = tri.Triangulation(mesh_coordinates[:, 0],
    #                                  mesh_coordinates[:, 1],
    #                                  triangles)

    # Plot the mesh
    # plt.triplot(triangulation)
    plt.xlabel("x")
    plt.ylabel("y")

    # Saving the mesh into the file
    fig = plt.gcf()
    fig.savefig(result_path + '-cells-%d-vertices-%d' % (n_cells, n_vert) + ".eps")
    plt.close('all')

def plot_decay_error_estimates_dofs(d, p, N, h, error, majorant, minorant, result_path):

    # Plot error and majorant results for p = 2
    plt.figure()
    plt.hold(True)
    plt.loglog(N, majorant, 'k-^', label="$maj$")
    plt.loglog(N, minorant, 'b-*', label="$[e]$")
    plt.loglog(N, error, 'r-s', label="$min$")

    # p -degree of the element
    # d - spacial dimension
    plt.loglog(N, 1 / (N**(p/d)), 'b:')

    plt.xlabel("h")
    plt.ylabel("rate")
    plt.legend()
    plt.grid(True)

    fig = plt.gcf()
    fig.savefig(result_path + "-maj-min-dofs.eps", format="eps")
    plt.close("all")

def plot_decay_error_dofs(d, p, N, h, e_h1, e_l2, e_linf, result_path):

    # Plot error and majorant results for p = 2
    plt.figure()
    plt.hold(True)
    plt.loglog(N, e_h1, 'k-^', label="$[e]$")
    plt.loglog(N, e_l2, 'b-*', label="$\|e\|_{L^2}$")
    plt.loglog(N, e_linf, 'r-s', label="$\|e\|_{L^{inf}}$")

    # p -degree of the element
    # d - spacial dimension
    plt.loglog(N, 1 / (N**(p / d)), 'b:', label="$r = p$")
    #plt.loglog(N, 1 / (N**((p - 1/2) / d)), 'g:', label="$r = p - 1/2$")

    plt.xlabel("h")
    plt.ylabel("rate")
    plt.legend()
    plt.grid(True)

    fig = plt.gcf()
    fig.savefig(result_path + "-l2-linf-error-dofs.eps", format="eps")
    plt.close("all")

def output_decay_error_estimates_summary(test_params, d, p, N, h, error, majorant, minorant):

    print '%------------------------------------------------------------------------------------%'
    print "% Rate of convergence summary"
    print '%------------------------------------------------------------------------------------%\n'

    print '       N &      [e] & r([e]) &      maj & r(maj) &      min & r(min) \\\\'
    print '-------------------------------------------------------------------------------'

    for i in range(0, len(h) - 1):
        h_i = 1 / (N[i] ** (1.0 / d))
        h_ip1 = 1 / (N[i + 1] ** (1.0 / d))

        maj_i = majorant[i]
        maj_ip1 = majorant[i + 1]

        min_i = minorant[i]
        min_ip1 = minorant[i + 1]

        e_i = error[i]
        e_ip1 = error[i + 1]

        rate_maj = ln(maj_ip1 / maj_i) / ln(h_ip1 / h_i)
        rate_min = ln(min_ip1 / min_i) / ln(h_ip1 / h_i)
        rate_e = ln(e_ip1 / min_i) / ln(h_ip1 / h_i)

        print ' %6d & %8.2e & %6.2f & %8.2e & %6.2f & %8.2e & %6.2f \\\\' \
              % (N[i], e_i, rate_e, maj_i, rate_maj, min_i, rate_min)
    print ""

    print "Chen paper: |ln(e) / ln(N)| = 1 indicates second convergence rate"
    print '        N &      [e] & ln([e])/ln(N) &      maj & r(maj) &      min & r(min) \\\\'
    print '-------------------------------------------------------------------------------'

    for i in range(0, len(h) - 1):
        maj_i = majorant[i]
        min_i = minorant[i]
        e_i = error[i]

        rate_maj = abs(ln(maj_i) / ln(N[i]))
        rate_min = abs(ln(min_i) / ln(N[i]))
        rate_e   = abs(ln(e_i)   / ln(N[i]))

        print ' %8.2e & %8.2e & %8.2f & %8.2e & %6.2f & %8.2e & %6.2f \\\\' \
              % (N[i], e_i, rate_e, maj_i, rate_maj, min_i, rate_min)
    print ""

    # Convergence rate in the adaptive refinement doesn't make sense wrt h
    if test_params["refinement_tag"] == 'uniform':

        print "% r = log2(maj_i / maj_ip1)"
        print ""
        print '        h &      [e] & r([e]) &      maj & r(maj) &      min & r(min) \\\\'
        print '-------------------------------------------------------------------------------'

        for i in range(0, len(h) - 1):

            maj_i = majorant[i]
            maj_ip1 = majorant[i + 1]

            min_i = minorant[i]
            min_ip1 = minorant[i + 1]

            e_i = error[i]
            e_ip1 = error[i + 1]

            # works !
            rate_maj = numpy.log2(maj_i / maj_ip1)
            rate_min = numpy.log2(min_i / min_ip1)
            rate_e = numpy.log2(e_i / e_ip1)

            print ' %8.2e & %8.2e & %6.2f & %8.2e & %6.2f & %8.2e & %6.2f \\\\' \
                  % (h[i], e_i, rate_e, maj_i, rate_maj, min_i, rate_min)
        print ""

        print "% r = ln(maj_i / maj_ip1) / ln(h_i / h_ip1) : "
        print ""
        print '        h &      [e] & r([e]) &      maj & r(maj) &      min & r(min) \\\\'
        print '-------------------------------------------------------------------------------'

        for i in range(0, len(h) - 1):
            h_i = h[i]
            h_ip1 = h[i + 1]

            maj_i = majorant[i]
            maj_ip1 = majorant[i + 1]

            min_i = minorant[i]
            min_ip1 = minorant[i + 1]

            e_i = error[i]
            e_ip1 = error[i + 1]

            rate_maj = ln(maj_i / maj_ip1) / ln(h_i / h_ip1)
            rate_min = ln(min_i / min_ip1) / ln(h_i / h_ip1)
            rate_e = ln(e_i / e_ip1) / ln(h_i / h_ip1)

            print ' %8.2e & %8.2e & %8.2f & %8.2e & %6.2f & %8.2e & %6.2f \\\\' \
                  % (h[i], e_i, rate_e, maj_i, rate_maj, min_i, rate_min)
        print ""

def output_decay_error_summary(test_params, d, p, N, h, e_array, e_l2_array, e_linf_array):

    print '%------------------------------------------------------------------------------------%'
    print "% Rate of convergence summary"
    print '%------------------------------------------------------------------------------------%\n'

    print ' N       &      [e] & r([e]) & \| e\|_{L2} & r(\| e\|_{L2}) & \| e\|_{Linf} & r(\| e\|_{Linf}) \\\\'
    print '-------------------------------------------------------------------------------'

    for i in range(0, len(h) - 1):
        h_i = 1 / (N[i] ** (1.0 / d))
        h_ip1 = 1 / (N[i + 1] ** (1.0 / d))

        e_l2_i = e_l2_array[i]
        e_l2_ip1 = e_l2_array[i + 1]

        e_linf_i = e_linf_array[i]
        e_linf_ip1 = e_linf_array[i + 1]

        e_energy_i = e_array[i]
        e_ip1 = e_array[i + 1]

        rate_e_l2 = ln(e_l2_ip1 / e_l2_i) / ln(h_ip1 / h_i)
        rate_e_linf = ln(e_linf_ip1 / e_linf_i) / ln(h_ip1 / h_i)
        rate_e = ln(e_ip1 / e_energy_i) / ln(h_ip1 / h_i)

        print ' %8d & %8.2e & %6.2f & %8.2e & %6.2f & %8.2e & %6.2f \\\\' \
              % (N[i], e_energy_i, rate_e, e_l2_i, rate_e_l2, e_linf_i, rate_e_linf)
    print ""

    print ""
    print "Chen paper: |ln(e) / ln(N)| = 1 indicates second convergence rate"
    print ' N       &      [e] & r([e]) & \| e\|_0 & r(\| e\|_0) & \| e\|_{inf} & r(\| e\|_{inf}) \\\\'
    print '-------------------------------------------------------------------------------'

    for i in range(0, len(h) - 1):
        e_l2_i = e_l2_array[i]
        e_linf_i = e_linf_array[i]
        e_energy_i = e_array[i]

        rate_e_l2 = abs(ln(e_l2_i) / ln(N[i]))
        rate_e_linf = abs(ln(e_linf_i) / ln(N[i]))
        rate_e = abs(ln(e_energy_i) / ln(N[i]))

        print ' %8d & %8.2e & %6.2f & %8.2e & %6.2f & %8.2e & %6.2f \\\\' \
              % (N[i], e_energy_i, rate_e, e_l2_i, rate_e_l2, e_linf_i, rate_e_linf)
    print ""

    # Convergence rate in the adaptive refinement doesn't make sense wrt h
    if test_params["refinement_tag"] == 'uniform':

        print "% r = log2(e_energy_i / e_ip1)"
        print ""
        print '        h &      [e] & r([e]) &    \|e\|_{L2} & r(\|e\|_{L2}) & \|e\|_{Linf} & r(\|e\|_{Linf}) \\\\'
        print '-------------------------------------------------------------------------------'

        for i in range(0, len(h) - 1):
            h_i = h[i]
            h_ip1 = h[i + 1]

            e_l2_i = e_l2_array[i]
            e_l2_ip1 = e_l2_array[i + 1]

            e_linf_i = e_linf_array[i]
            e_linf_ip1 = e_linf_array[i + 1]

            e_energy_i = e_array[i]
            e_ip1 = e_array[i + 1]

            # works !
            rate_e_l2 = numpy.log2(e_l2_i / e_l2_ip1)
            rate_e_linf = numpy.log2(e_linf_i / e_linf_ip1)
            rate_e = numpy.log2(e_energy_i / e_ip1)

            print ' %8.2e & %8.2e & %6.2f & %8.2e & %6.2f & %8.2e & %6.2f \\\\' \
                  % (h_i, e_energy_i, rate_e, e_l2_i, rate_e_l2, e_linf_i, rate_e_linf)
        print ""

        print "% r = ln(e_energy_i / e_ip1) / ln(h_i / h_ip1) : "
        print ""
        print '        h &      [e] & r([e]) &     \|e\|_{L2} & r(\|e\|_{L2}) & \|e\|_{Linf} & r(\|e\|_{Linf}) \\\\'
        print '-------------------------------------------------------------------------------'

        for i in range(0, len(h) - 1):
            h_i = h[i]
            h_ip1 = h[i + 1]

            e_l2_i = e_l2_array[i]
            e_l2_ip1 = e_l2_array[i + 1]

            e_linf_i = e_linf_array[i]
            e_linf_ip1 = e_linf_array[i + 1]

            e_energy_i = e_array[i]
            e_ip1 = e_array[i + 1]

            rate_e_l2 = ln(e_l2_ip1 / e_l2_i) / ln(h_ip1 / h_i)
            rate_e_linf = ln(e_linf_ip1 / e_linf_i) / ln(h_ip1 / h_i)
            rate_e = ln(e_ip1 / e_energy_i) / ln(h_ip1 / h_i)

            print ' %8.2e & %8.2e & %6.2f & %8.2e & %6.2f & %8.2e & %6.2f \\\\' \
                  % (h_i, e_energy_i, rate_e, e_l2_i, rate_e_l2, e_linf_i, rate_e_linf)
        print ""

def output_error_estimates_summary(N, e_maj, maj, e_min, min):

    print '%------------------------------------------------------------------------------------%'
    print "% Efficiency indices summary"
    print '%------------------------------------------------------------------------------------%\n'

    print ' REF \# &     DOFs & [e]_{maj} &      maj & i_eff(maj) & [e]_{min} &      min & i_eff(min) \\\\'
    print '---------------------------------------------------------------------------------'
    for i in range(0, len(N)):
        print ' %6d & %8d & %9.2e & %8.2e & %10.2f & %9.2e & %8.2e & %10.2f \\\\' \
              % (i, N[i], e_maj[i], maj[i], maj[i]/e_maj[i], e_min[i], min[i], min[i]/e_min[i])
    print ""

def plot_decay_majorant_dofs_v_of_deg(d, p, N, maj, result_path):

    plt.figure()
    plt.hold(True)
    plt.loglog(N, maj, 'g-*')

    # p -degree of the element
    # d - spacial dimension
    plt.loglog(N, 1 / (N ** (p / d)), 'b:')


    plt.xlabel("h")
    plt.ylabel("M^2 P%d" % (p))
    plt.grid(True)

    fig = plt.gcf()
    fig.savefig(result_path + "-dofs.eps", format="eps")
    plt.close('all')

def save_decay_to_mat_file(dofs, error, majorant, result_path):
    scipy.io.savemat(result_path, {'dofs': dofs, 'error': error, 'majorant': majorant})


def plot_mesh_3d(mesh, result_path):

    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()

    d = mesh.geometry().dim()
    mesh_coordinates = mesh.coordinates()

    # "-num-cells-%d-num-vertices-%d.eps"
    x = mesh_coordinates[:, 0]
    y = mesh_coordinates[:, 1]
    z = mesh_coordinates[:, 2]
    scipy.io.savemat(result_path + '-cells-%d-vertices-%d-xyz-triangulation' % (num_cells, num_vert),
                     {'x': x,
                      'y': y,
                      'z': z})
    plt.close('all')

def plot_carpet_2d(mesh, f, result_path):
    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()

    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates().reshape((num_vert, d))
    triangles = numpy.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)

    # Get the z values as face colors for each triangle(midpoint)
    fig = plt.figure()

    # Get the z values for each vertex
    plt.figure()
    plt.xlabel('x')
    plt.ylabel('y')
    z = numpy.asarray([f(point) for point in mesh_coordinates])
    plt.tripcolor(triangulation, z, cmap=cm.coolwarm, edgecolors='k')
    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))
    plt.close('all')

# Plotting the scalar function
def plot_function(f, mesh, dim, result_path):

    if dim == 2:
        plot_function_2d(mesh, f, result_path)
        plot_function_3d(mesh, f, result_path)

def plot_function_1d(mesh, f, result_path):

    num_cells = mesh.num_cells()
    num_vert = mesh.num_vertices()

    plotting.plot(f, mesh=mesh)
    plt.xlabel('x')
    plt.ylabel('u')
    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))
    plt.close('all')

    """
    num_cells = mesh.num_cells()
    num_vert = mesh.num_vertices()

    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates()
    z = np.asarray([f(point) for point in mesh_coordinates])

    plt.figure()
    plt.plot(mesh_coordinates, z, 'b.-')
    plt.xlabel('x')
    plt.ylabel('u')

    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))
    plt.close('all')
    """

def plot_function_2d(mesh, f, result_path):

    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()

    d = mesh.geometry().dim()

    plotting.plot(f, cmap=cm.BrBG, edgecolors='k', mesh=mesh)

    # Create the triangulation
    #mesh_coordinates = mesh.coordinates().reshape((num_vert, d))
    #triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    #triangulation = tri.Triangulation(mesh_coordinates[:, 0],
    #                                  mesh_coordinates[:, 1],
    #                                  triangles)

    '''
    # Get the z values as face colors for each triangle(midpoint)
    plt.figure()
    zfaces = np.asarray([f(cell.midpoint()) for cell in cells(mesh)])
    plt.tripcolor(triangulation, facecolors=zfaces, edgecolors='k')
    fig = plt.gcf()

    fig.savefig(result_path + '-zfaces-cells-%d-vertices-%d.eps' %(num_cells, num_vert))
    '''

    fig = plt.gcf()
    fig.savefig(result_path + '-z-cells-%d-vertices-%d.eps' % (num_cells, num_vert))
    plt.close('all')

def plot_function_3d(mesh, f, result_path):

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm

    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()

    d = mesh.geometry().dim()

    mesh_coordinates = mesh.coordinates().reshape((num_vert, d))
    triangles = numpy.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)
    z = numpy.asarray([f(point) for point in mesh_coordinates])

    # change the font
    matplotlib.rcParams.update({'font.size': 10,
                                'font.family': 'serif'})

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    '''
    # customize ticks on the axises
    start, end = ax.get_zlim()
    ax.zaxis.set_ticks(numpy.arange(start, end, (end - start) / 3))
    ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(numpy.arange(start, end, (end - start) / 3))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(numpy.arange(start, end, (end - start) / 3))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))
    '''

    # add labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u')

    # ax.plot_trisurf(triangulation, z, cmap=cm.gray, linewidth=0.2) RdGy
    ax.plot_trisurf(triangulation, z, cmap=cm.BrBG, linewidth=0.2)
    # ax.plot_trisurf(triangulation, z, cmap=cm.jet, linewidth=0.2)
    fig.savefig(result_path + '-cells-%d-vertices-%d.eps' % (num_cells, num_vert))
    plt.close('all')

