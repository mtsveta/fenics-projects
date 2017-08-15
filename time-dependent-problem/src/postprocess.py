from dolfin.cpp.mesh import cells

__author__ = 'svetlana'

import scipy.io
import os, sys
import numpy
from dolfin import *

import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.ticker as ticker
import matplotlib
from matplotlib import cm
#from mpl_toolkits.mplot3d import axes3d, Axes3D
#from itertools import product, combinations
import plotting

import plotslopes


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


def output_problem_characteristics(test_num, u_expr, grad_u_x, f_expr,
                                   t_T, domain, dim,
                                   nt, num_cells, num_vertices,
                                   v_approx_deg, flux_approx_deg, delta):
    print "test: ", test_num
    print "domain: ", domain
    print "u = ", u_expr
    print "f = ", f_expr

    print 'mesh parameters: time steps = %d, num_cells = %d, num_vertices = %d' % (nt, num_cells, num_vertices)
    print "T = ", t_T
    print "space func approx_deg = ", v_approx_deg
    print "flux approx_deg = ", flux_approx_deg
    print "delta = ", delta


def output_time_layer_result_error_and_majorant(t_k, t_kp1, error_k, majorant_k, i_eff_maj_k):
    print 'error     = %8.2e ' % (error_k)
    print 'majorant  = %8.2e' % (majorant_k)
    print 'i_eff_maj = %.4f' % (i_eff_maj_k)

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

    #plotting.plot(mesh)
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


def plot_mesh_3d(mesh, result_path):
    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()

    plt.figure()
    plotting.plot(mesh, cmap=cm.coolwarm)
    plt.xlabel("x")
    plt.ylabel("y")

    # Saving the mesh into the file
    fig = plt.gcf()
    fig.savefig(result_path + '-cells-%d-vertices-%d' % (num_cells, num_vert) + ".eps")
    plt.close('all')

    d = mesh.geometry().dim()
    mesh_coordinates = mesh.coordinates()

    # "-num-cells-%d-num-vertices-%d.eps"
    x = mesh_coordinates[:, 0]
    y = mesh_coordinates[:, 1]
    z = mesh_coordinates[:, 2]
    scipy.io.savemat(result_path + '-cells-%d-vertices-%d-xyz-triangulation.eps' % (num_cells, num_vert),
                     {'x': x,
                      'y': y,
                      'z': z})

def save_results_to_mat_file(error_distr, maj_distr,
                             error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_max_k, h_min_k,
                             result_path):
    scipy.io.savemat(result_path,
                     {'error_distr': error_distr, 'maj_distr': maj_distr,
                      'error_k': error_k, 'majorant_k': majorant_k,
                      'i_eff_maj_k': i_eff_maj_k, 'h_max_k': h_max_k, 'h_min_k': h_min_k,
                      'error_d_k': error_d_k, 'majorant_d_k': majorant_d_k})


def save_distr_to_mat_file(error_distr, maj_distr, result_path):
    scipy.io.savemat(result_path,
                     {'error_distr': error_distr, 'maj_distr': maj_distr})


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
    n = mesh.num_vertices()
    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates().reshape((n, d))
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)

    # Get the z values as face colors for each triangle(midpoint)
    plt.figure()
    zfaces = np.asarray([f(cell.midpoint()) for cell in cells(mesh)])
    plt.tripcolor(triangulation, facecolors=zfaces, edgecolors='k')
    fig = plt.gcf()
    fig.savefig(result_path + "-faces.eps")

    # Get the z values for each vertex
    plt.figure()
    z = np.asarray([f(point) for point in mesh_coordinates])
    plt.tripcolor(triangulation, z, edgecolors='k')
    plt.xlabel('x')
    plt.ylabel('u')

    fig = plt.gcf()
    fig.savefig(result_path + "-vertices.eps")


def plot_function_3d(mesh, f, result_path):
    from matplotlib import cm

    n_vert = mesh.num_vertices()
    n_cells = mesh.num_cells()

    d = mesh.geometry().dim()

    mesh_coordinates = mesh.coordinates().reshape((n_vert, d))
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)
    z = np.asarray([f(point) for point in mesh_coordinates])

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot_trisurf(triangulation, z, cmap=cm.coolwarm, linewidth=0.2)
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" %(n_cells, n_vert))


def plot_function(f, mesh, dim, result_path):
    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()

    plot(f, interactive=True)
    if dim == 1:
        plot_function_1d(mesh, num_cells, num_vert, f, result_path)
    elif dim == 2:
        plot_function_2d(mesh, f, result_path)
        plot_function_3d(mesh, f, result_path)


def document_results_without_refinement(mesh, e_distr, m_distr, error_k, majorant_k, i_eff_maj_k, error_d_k,
                                        majorant_d_k, h_max_k, h_min_k, t_T, k, nt, project_path, results_folder,
                                        results_info):
    if k == nt - 1:
        plot_integral_error_majorant(numpy.linspace(0, t_T, k), error_k[0:k], majorant_k[0:k],
                                     project_path + results_folder + 'integral-e-maj-no-refinement' + results_info)
        # Save the results in mat-file
        save_results_to_mat_file(e_distr, m_distr, error_k, majorant_k, i_eff_maj_k, error_d_k, majorant_d_k, h_max_k,
                                 h_min_k,
                                 project_path + results_folder + 'no-refinement-final-results' + results_info)
        plot_histogram(mesh, e_distr, m_distr,
                       project_path + results_folder + 'e-maj-distr-hist' + results_info)


def document_results_with_refinement(problem_params, mesh, e_distr, m_distr, error_k, majorant_k, i_eff_maj_k,
                                     error_d_k, majorant_d_k, h_max_k, h_min_k, t_T, k, nt, refinement_tag,
                                     project_path, results_folder, results_info, dim):
    if k + 1 == nt:
        plot_integral_error_majorant(numpy.linspace(0, t_T, k), error_k[0:k], majorant_k[0:k],
                                     project_path + results_folder + 'int-e-maj-' + refinement_tag + results_info)
        # Save the results in mat-file
        save_results_to_mat_file(e_distr, m_distr, error_k, majorant_k, i_eff_maj_k, h_max_k, h_min_k, error_d_k,
                                 majorant_d_k,
                                 project_path + results_folder + 'final-results-' + refinement_tag + results_info)
    results_path = project_path + results_folder + 'mesh-' + refinement_tag + results_info
    if dim == 2:    plot_mesh(mesh, results_path)
    elif problem_params['refinement_tag'] == "adaptive" and dim == 3: plot_mesh_3d(mesh, results_path)

def plot_function_1d(mesh, num_cells, num_verts, f, result_path):

    # Create the triangulation
    mesh_coordinates = mesh.coordinates()
    z = numpy.asarray([f(point) for point in mesh_coordinates])

    plt.figure()
    plt.plot(mesh_coordinates, z, 'b.-')
    plt.xlabel('x')
    plt.ylabel('u')

    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_verts))


def plot_carpet_2d(mesh, f, result_path):

    num_cells = mesh.num_cells()
    num_verts = mesh.num_vertices()
    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates().reshape((num_verts, d))
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
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_verts))

def plot_function_2d(mesh, num_cells, num_verts, f, result_path):
    n = mesh.num_vertices()
    d = mesh.geometry().dim()

    plotting.plot(f, cmap=cm.BrBG, edgecolors='k', mesh=mesh)

    fig = plt.gcf()
    fig.savefig(result_path + "-faces.eps")

    fig = plt.gcf()
    fig.savefig(result_path + '-z-cells-%d-vertices-%d.eps' % (num_cells, num_verts))
    plt.close('all')

def plot_function_3d(mesh, f, result_path):
    from mpl_toolkits.mplot3d import Axes3D


    n_vert = mesh.num_vertices()
    n_cell = mesh.num_cells()

    d = mesh.geometry().dim()

    mesh_coordinates = mesh.coordinates().reshape((n_vert, d))
    triangles = numpy.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)
    z = numpy.asarray([f(point) for point in mesh_coordinates])

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(triangulation, z, cmap=cm.coolwarm, linewidth=0.2)

    # change the font
    matplotlib.rcParams.update({'font.size': 16,
                                'font.family': 'serif'})

    # customize ticks on the axises
    start, end = ax.get_zlim()
    ax.zaxis.set_ticks(numpy.arange(start, end, (end - start) / 5))
    ax.zaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

    start, end = ax.get_xlim()
    ax.xaxis.set_ticks(numpy.arange(start, end, (end - start) / 5))
    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

    start, end = ax.get_ylim()
    ax.yaxis.set_ticks(numpy.arange(start, end, (end - start) / 5))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

    # add labels
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('u')

    fig.savefig(result_path + '-cells-%d-vertices-%d.eps' %(n_cell, n_vert))

def plot_solution(f, mesh, dim, result_path):
    if dim == 2:
        plot_function_3d(mesh, f, result_path)

def construct_result_folder_name(dim, test_num, problem_params):
    results_folder_name = '/out/%dd/test-%d/v-P%d-y-RT%d/' \
                          % (dim, test_num,
                             problem_params["v_approx_order"],
                             problem_params["flux_approx_order"])
    if problem_params["solving_strategy_tag"] == "without-refinement":
        results_folder_name += 'no-ref/'
    else:
        if problem_params["refinement_tag"] == 'uniform':
            results_folder_name += 'uniform-ref/'
        else:
            results_folder_name += problem_params["marking_tag"] + '/'
    return results_folder_name


