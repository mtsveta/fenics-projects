__author__ = 'svetlana'

import scipy.io
import os
import numpy
from dolfin import *
from matplotlib import cm
from dolfin.cpp.mesh import cells

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import plotslopes
import matplotlib
import plotting

def allocate_array(n):
    return numpy.array([0.0 for i in range(n)])

def allocate_array_2d(m, n):
    return numpy.array([[0.0 for j in range(n)] for i in range(m)])

def save_distr_to_mat_file(eta_omega_cell, E_DWR_cell, J_e_cell, M_cell, result_path):
    scipy.io.savemat(os.path.split(os.path.abspath(os.path.dirname(__file__)))[0] + result_path,
                     {'eta_omega_cell': eta_omega_cell,
                      'E_DWR_cell': E_DWR_cell,
                      'J_e_cell': J_e_cell,
                      'M_cell': M_cell})

def save_marking_to_mat_file(test_num, num_cell, V, marking, path):
    scipy.io.savemat(os.path.split(os.path.abspath(os.path.dirname(__file__)))[0] + path,
                     {'marking': marking})

def output_problem_characteristics(test_num, u0_expr, uD_expr, f_expr, A_expr, lmbd_min,
                                   t_T, domain, dim,
                                   nt, num_cells, num_vertices,
                                   v_approx_deg, flux_approx_deg, delta):

    print "test: ", test_num
    print "domain: ", domain
    print "u0 = ", u0_expr
    print "uD = ", uD_expr
    print "f = ", f_expr
    print "A = ", A_expr
    print "lmbd_min = ", lmbd_min

    print 'mesh parameters: time steps = %d, num_cells = %d, num_vertices = %d' % (nt, num_cells, num_vertices)
    print "T = ", t_T
    print "space func approx_deg = ", v_approx_deg
    print "flux approx_deg = ", flux_approx_deg
    print "delta = ", delta



def output_time_layer_result_error_and_majorant(t_k, t_kp1, error_k, majorant_k, i_eff_maj_k):
    print 'error     = %8.2e ' % (error_k)
    print 'majorant  = %8.2e' % (majorant_k)
    print 'i_eff_maj = %.4f' % (i_eff_maj_k)

def log_results_summary(err, maj, i_eff):
    print '%------------------------------------------------------------------------------------%'
    print "% Efficiency indices summary"
    print '%------------------------------------------------------------------------------------%\n'

    print ' k &     [e]^2 &      maj^2 & i_eff(maj)\\\\'
    print '---------------------------------------------------------------------------------'
    for i in range(0, len(err)):
        print ' %6d & %9.2e & %8.2e & %10.2f \\\\' \
              % (i, err[i], maj[i], i_eff[i])
    print ""

def construct_result_tag(test_num, nx0, nx1, nx2, nt, k):
    file_info = '-test-%d-nt-%d-%d-%d-%d-time-%d' \
                % (test_num, nt, nx0, nx1, nx2, k)
    return file_info

def plot_histogram(mesh, e_distr, m_distr, result_path):

    num_cells = mesh.num_cells()
    num_vert = mesh.num_vertices()
    cells = numpy.linspace(0, num_cells, num_cells)

    plt.figure()
    plt.bar(cells, e_distr, align="center", facecolor="red")
    plt.plot(cells, m_distr, 'g.-')
    plt.xlabel('cells')
    plt.ylabel('e distr, m distr')

    fig = plt.gcf()

    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))
    plt.close('all')

def plot_mesh(mesh, result_path):

    plt.figure()
    plt.hold(True)

    # change the font
    matplotlib.rcParams.update({'font.size': 16,
                                'font.family': 'serif'})


    n_vert = mesh.num_vertices()
    n_cells = mesh.num_cells()

    d = mesh.geometry().dim()

    # Create the triangulation
    plotting.plot(mesh)

    # Plot the mesh
    plt.xlabel("x")
    plt.ylabel("y")

    fig = plt.gcf()
    fig.savefig(result_path + '-cells-%d-vertices-%d' % (n_cells, n_vert) + ".eps")

def save_results_to_mat_file(error_distr, maj_distr,
                             error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_k,
                             result_path):
    scipy.io.savemat(result_path,
                     {'error_distr': error_distr, 'maj_distr': maj_distr,
                      'error_k': error_k, 'majorant_k': majorant_k,
                      'i_eff_maj_k': i_eff_maj_k, 'h_k': h_k,
                      'error_d_k': error_d_k, 'majorant_d_k': majorant_d_k})
def save_distr_to_mat_file(error_distr, maj_distr, maj_f_distr, result_path):
    scipy.io.savemat(result_path,
                     {'error_distr': error_distr, 'maj_distr': maj_distr, 'maj_f_distr': maj_f_distr})

def plot_integral_error_majorant(time, error_k, majorant_k, result_path):
    plt.figure()
    plt.plot(time, error_k, 'y.-')
    plt.plot(time, majorant_k, 'b.-')
    plt.xlabel('t')
    plt.ylabel('e, maj')
    #plt.show()

    fig = plt.gcf()
    fig.savefig(result_path + ".eps")

def get_project_path():
    (project_path, src_folder) = os.path.split(os.path.abspath(os.path.dirname(__file__)))
    #(project_path, src_folder) = os.path.split(os.path.abspath(src_path))
    return project_path

def create_results_folder(directory):
    full_directory = get_project_path() + directory

    if not os.path.exists(full_directory):
        os.makedirs(full_directory)

    #sys.stdout = open(full_directory + "log.dat", "w+")

def plot_function_1d(mesh, f, result_path):

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


def plot_decay_error_majorant(deg, h, error, majorant, result_path):

    # Plot error and majorant results for deg = 2
    plt.figure()
    plt.hold(True)
    plt.loglog(h, majorant, 'k-^')
    plt.loglog(h, error, 'k-s')
    if deg == 1:
        plotslopes.slope_marker((h[1], error[1]), (2, 1))
    elif deg == 2:
        plotslopes.slope_marker((h[1], error[1]), (3, 1))

    plt.legend(["M^2", "|| e ||^2"], "lower right")
    plt.xlabel("h")
    plt.ylabel("|| e ||^2, M^2")
    plt.grid(True)

    fig = plt.gcf()
    fig.savefig(result_path + ".eps")

def plot_decay_majorant_v_of_deg(h, maj, deg, result_path):

    plt.figure()
    plt.hold(True)
    plt.loglog(h, maj, 'b-*')
    plotslopes.slope_marker((h[1], maj[1]), (1 + deg, 1))
    plt.legend(["M^2"], "lower right")
    plt.xlabel("h")
    plt.ylabel("M^2 P%d" % (deg))
    plt.grid(True)

    fig = plt.gcf()
    fig.savefig(result_path + ".eps")

def save_decay_to_mat_file(h, error, majorant, result_path):
    scipy.io.savemat(result_path, {'h': h, 'error': error, 'majorant': majorant})

def plot_function_2d(mesh, f, result_path):

    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()
    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates().reshape((num_vert, d))
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)

    # Get the z values as face colors for each triangle(midpoint)
    plt.figure()
    zfaces = np.asarray([f(cell.midpoint()) for cell in cells(mesh)])
    plt.tripcolor(triangulation, facecolors=zfaces, edgecolors='k')
    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))


    # Get the z values for each vertex
    plt.figure()
    z = np.asarray([f(point) for point in mesh_coordinates])
    plt.tripcolor(triangulation, z, edgecolors='k')
    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))


def plot_function_3d(mesh, f, result_path):

    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()

    d = mesh.geometry().dim()

    mesh_coordinates = mesh.coordinates().reshape((num_vert, d))
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)
    z = np.asarray([f(point) for point in mesh_coordinates])


    fig = plt.figure()
    ax = fig.gca(projection='3d')

    #ax.plot_trisurf(triangulation, z, cmap=cm.jet, linewidth=0.2)
    ax.plot_trisurf(triangulation, z, cmap=cm.coolwarm, linewidth=0.2)
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))


def plot_function(f, mesh, dim, result_path):

    plot(f, interactive=True)
    if dim == 1:
        plot_function_1d(mesh, f, result_path)
    elif dim == 2:
        plot_function_2d(mesh, f, result_path)
        plot_function_3d(mesh, f, result_path)


def document_results_without_refinement(mesh, e_distr, m_distr, error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_k, t_T, k, nt, project_path, results_folder, results_info):
    if k == nt-1:
        plot_integral_error_majorant(numpy.linspace(0, t_T, k), error_k[0:k], majorant_k[0:k],
                                                     project_path + results_folder + 'integral-e-maj-no-refinement' + results_info)
        # Save the results in mat-file
        save_results_to_mat_file(e_distr, m_distr, error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_k,
                                             project_path + results_folder + 'no-refinement-final-results-' + results_info)
        plot_histogram(mesh, e_distr, m_distr,
                       project_path + results_folder + 'e-maj-distr-hist' + results_info)
    plot_mesh(mesh, project_path + results_folder + 'mesh-' + results_info)

def document_results_with_refinement(problem_params, mesh, e_distr, m_distr, maj_distr,
                                     error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k,
                                     h_k, t_T, k, nt, refinement_tag, project_path, results_folder, results_info, dim):
    if k+1 == nt:
        plot_integral_error_majorant(numpy.linspace(0, t_T, k), error_k[0:k], majorant_k[0:k],
                                                 project_path + results_folder + 'int-e-maj-' + refinement_tag + results_info)
        # Save the results in mat-file
        save_results_to_mat_file(e_distr, m_distr, error_k, majorant_k, i_eff_maj_k, h_k, error_d_k, majorant_d_k,
                                             project_path + results_folder + 'final-results-' + refinement_tag + results_info)

    if (problem_params['refinement_tag'] == "uniform"  or problem_params['refinement_tag'] == "adaptive") and dim == 2:
        plot_mesh(mesh, project_path + results_folder + 'mesh-' + refinement_tag + results_info)

    if problem_params['solving_strategy_tag'] == "with-refinement":
        # postprocess.plot_histogram(mesh, e_distr, m_distr,
        #                           project_path + results_folder + 'e-maj-distr-hist-' + problem_params['refinement_criteria_tag'] + '-' + problem_params['marking_tag'] + '-marking-' + results_info)
        save_distr_to_mat_file(e_distr, m_distr, maj_distr,
                               project_path + results_folder + 'e-maj-distr-hist-' + refinement_tag + results_info)

def plot_function_1d(mesh, f, result_path):
    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()

    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates()
    z = numpy.asarray([f(point) for point in mesh_coordinates])

    plt.figure()
    plt.plot(mesh_coordinates, z, 'b.-')
    plt.xlabel('x')
    plt.ylabel('u')

    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps")


def plot_function_2d(mesh, f, result_path):

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
    plt.figure()
    zfaces = numpy.asarray([f(cell.midpoint()) for cell in cells(mesh)])
    plt.tripcolor(triangulation, facecolors=zfaces, edgecolors='k')
    fig = plt.gcf()

    fig.savefig(result_path + '-cells-%d-vertices-%d.eps' %(num_cells, num_vert))

    # Get the z values for each vertex
    plt.figure()
    z = numpy.asarray([f(point) for point in mesh_coordinates])
    plt.tripcolor(triangulation, z, edgecolors='k')
    fig = plt.gcf()

    fig.savefig(result_path + '-cells-%d-vertices-%d.eps' %(num_cells, num_vert))

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
    matplotlib.rcParams.update({'font.size': 16,
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

    #ax.plot_trisurf(triangulation, z, cmap=cm.gray, linewidth=0.2) RdGy
    ax.plot_trisurf(triangulation, z, cmap=cm.coolwarm, linewidth=0.2)
    #ax.plot_trisurf(triangulation, z, cmap=cm.jet, linewidth=0.2)
    fig.savefig(result_path + '-cells-%d-vertices-%d.eps' %(num_cells, num_vert))


def plot_solution(f, mesh, dim, result_path):

    if dim == 2:
        #plot_function_2d(mesh, f, result_path)
        plot_function_3d(mesh, f, result_path)

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

    fig = plt.figure()

    # add labels
    plt.figure()
    plt.xlabel('x')
    plt.ylabel('y')

    # change the font
    matplotlib.rcParams.update({'font.size': 16,
                                'font.family': 'serif'})


    z = numpy.asarray([f(point) for point in mesh_coordinates])
    plt.tripcolor(triangulation, z, cmap=cm.coolwarm, edgecolors='k')
    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))
