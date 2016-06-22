from dolfin.cpp.mesh import cells

__author__ = 'svetlana'

import scipy.io
import os, sys
import numpy

import matplotlib.pylab as fh
import matplotlib.pylab as plt
import matplotlib.pyplot as plt
import matplotlib.tri as tri

from dolfin import *
import plotslopes
import math
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib.ticker as ticker
import matplotlib
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d, Axes3D
from itertools import product, combinations

def construct_result_tag(test_num, ref_num, nx0, nx1, nx2, num_cells, num_vertices):
    file_info = '-test-num-%d-ref-inter-%d-%d-%d-%d-num-cells-%d-num-vertices-%d' \
                % (test_num, ref_num, nx0, nx1, nx2, num_cells, num_vertices)
    return file_info

def allocate_array(n):
    return numpy.array([0.0 for i in range(n)])

def allocate_array_2d(m, n):
    return numpy.array([[0.0 for j in range(n)] for i in range(m)])

def plot_histogram(mesh, e_distr, m_distr, result_path):

    num_cells = mesh.num_cells()
    cells = numpy.linspace(0, num_cells, num_cells)

    plt.figure()
    plt.plot(cells, e_distr, 'r.-')
    plt.plot(cells, m_distr, 'g.-')
    plt.xlabel('cells')
    plt.ylabel('e distr, m distr')
    #plt.show()

    fig = plt.gcf()
    fig.savefig(result_path + '.eps', format="eps")


def plot_integral_error_majorant(time, error_k, majorant_k, result_path):
    plt.figure()
    plt.plot(time, error_k, 'y.-')
    plt.plot(time, majorant_k, 'b.-')
    plt.xlabel('t')
    plt.ylabel('e, maj')
    #plt.show()

    fig = plt.gcf()
    fig.savefig(result_path+ ".eps", format="eps")

def get_project_path():
    (project_path, src_folder) = os.path.split(os.path.abspath(os.path.dirname(__file__)))
    #(project_path, src_folder) = os.path.split(os.path.abspath(src_path))
    return project_path

def create_results_folder(directory):
    full_directory = get_project_path() + directory

    if not os.path.exists(full_directory):
        os.makedirs(full_directory)
    #sys.stdout = open(full_directory + "log.dat", "w+")


def plot_mesh(mesh, result_path):

    plt.figure()
    plt.hold(True)

    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()
    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates().reshape((num_vert, d))
    triangles = numpy.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)

    # Plot the mesh

    # change the font
    matplotlib.rcParams.update({'font.size': 16,
                                'font.family': 'serif'})

    plt.triplot(triangulation)
    plt.xlabel("x")
    plt.ylabel("y")

    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))

def save_results_to_mat_file(error_distr, maj_distr, error, majorant, i_eff_maj, result_path):
    scipy.io.savemat(result_path,
                     {'error_distr': error_distr, 'maj_distr': maj_distr,
                      'error': error, 'majorant': majorant, 'i_eff_maj': i_eff_maj})


def plot_refinement_error_majorant(h, error, majorant, result_path):

    plt.figure()
    plt.plot(numpy.log(h), numpy.log(error), 'y.-')
    plt.plot(numpy.log(h), numpy.log(majorant), 'b.-')
    plt.xlabel('log h')
    plt.ylabel('log e, log maj')
    #plt.show()

    fig = plt.gcf()
    fig.savefig(result_path + ".eps", format="eps")

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
    fig.savefig(result_path + ".eps", format="eps")

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
    fig.savefig(result_path + ".eps", format="eps")

def save_decay_to_mat_file(h, error, majorant, result_path):
    scipy.io.savemat(result_path, {'h': h, 'error': error, 'majorant': majorant})


def plot_mesh_3d(mesh, result_path):
    #plt.figure()
    #plt.hold(True)

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

def plot_function_3d(mesh, f, result_path):

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm

    n = mesh.num_vertices()
    d = mesh.geometry().dim()

    mesh_coordinates = mesh.coordinates().reshape((n, d))
    triangles = numpy.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)
    z = numpy.asarray([f(point) for point in mesh_coordinates])


    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ax.plot_trisurf(triangulation, z, cmap=cm.jet, linewidth=0.2)
    fig.savefig(result_path + "-3d.eps", format="eps")

def plot_function_2d(mesh, f, result_path):
    n = mesh.num_vertices()
    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates().reshape((n, d))
    triangles = numpy.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)
    '''
    # Get the z values as face colors for each triangle(midpoint)
    fig = plt.figure()
    zfaces = numpy.asarray([f(cell.midpoint()) for cell in cells(mesh)])
    plt.tripcolor(triangulation, facecolors=zfaces, edgecolors='k')
    plt.xlabel('x', fontsize=8)
    plt.ylabel('y', fontsize=8)

    fig = plt.gcf()
    fig.savefig(result_path + "-faces.eps")
    '''
    # Get the z values for each vertex
    plt.figure()
    z = numpy.asarray([f(point) for point in mesh_coordinates])
    plt.tripcolor(triangulation, z, edgecolors='k')
    fig = plt.gcf()
    fig.savefig(result_path + "-vertices.eps")

def plot_function(f, mesh, dim, result_path):

    if dim == 2:
        plot_function_2d(mesh, f, result_path)
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

