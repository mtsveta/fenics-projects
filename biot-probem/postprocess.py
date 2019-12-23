__author__ = 'svetlana'

from dolfin import *

import scipy.io
import os, sys
import numpy as np
import math

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import matplotlib
import plotting

# -------------------------------------------------------------------------------------------------------------------- #
# IO functions
# -------------------------------------------------------------------------------------------------------------------- #
# Function returning the path of the project
def get_project_path():
    (project_path, src_folder) = os.path.split(os.path.abspath(os.path.dirname(__file__)))
    #(project_path, src_folder) = os.path.split(os.path.abspath(src_path))
    return project_path

# Constructing the name of the result folder based on the problem parameters
def construct_result_folder_name(dim, test_num, problem_params):
    if problem_params["coupling_approach"] == 'iterative':
        results_folder_name = '/out-3d/%dd/test-%d/iterative/test-%d-v-P%d-p-P%d-flux-RT%d-tau-P%d' \
                              '-nt' \
                              '-%d' \
                              '-nit-%d-meshes-%d/' \
                              % (dim, test_num, test_num,
                                 problem_params["u_approx_order"],
                                 problem_params["p_approx_order"],
                                 problem_params["flux_approx_order"],
                                 problem_params["stress_approx_order"],
                                 problem_params["time_steps"],
                                 problem_params["iter_num"],
                                 problem_params["mesh_resolution"])
    else:
        results_folder_name = '/out/%dd/test-%d/full-implicit/test-%d-v-P%d-p-P%d-flux-RT%d-tau-P%d-nt-%d-meshes-%d/' \
                              % (dim, test_num, test_num,
                                 problem_params["u_approx_order"],
                                 problem_params["p_approx_order"],
                                 problem_params["flux_approx_order"],
                                 problem_params["stress_approx_order"],
                                 problem_params["time_steps"],
                                 problem_params["mesh_resolution"])
    '''
    if problem_params["solving_strategy_tag"] == "without-refinement":
        results_folder_name += 'no-ref/'
    else:
        if problem_params["refinement_tag"] == 'uniform':
            results_folder_name += 'uniform-ref/'
        else:
            results_folder_name += problem_params["marking_tag"] + '/'
    '''
    return results_folder_name

# Function creating the result folder
def create_results_folder(test_params, directory):
    # Creat the full directory of the result folder
    full_directory = get_project_path() + directory

    # If derictory doesn't exist
    if not os.path.exists(full_directory):
        os.makedirs(full_directory)
    # Write into the logging path
    if (test_params["output"]=="file"):
        sys.stdout = open(full_directory + "log-results.txt", "w+")

# Saving data into the mat-file
def save_to_mat_file(data, data_tag, result_path):
    scipy.io.savemat(os.path.split(os.path.abspath(os.path.dirname(__file__)))[0] + result_path,
                     {data_tag: data})
# Plotting the problem data
def output_problem_characteristics(test_num, mesh, problem_data, domain_params, material_params, test_params, func_spaces):

    print("test: ", test_num)
    print("domain: ", domain_params['domain_type'])
    print("l_x: ", domain_params['l_x'])
    print("l_y: ", domain_params['l_y'])
    print("l_z: ", domain_params['l_z'])
    print("gdim: ", domain_params['gdim'])

    print("u_expr: ", problem_data['u_expr'])
    print("p_expr: ", problem_data['p_expr'])
    print("t_T: ", problem_data['t_T'])

    print("alpha: ", material_params['alpha'])
    print("beta: ", material_params['beta'])

    print("lmbda: ", material_params['lmbda'])
    print("mu: ", material_params['mu'])
    print("kappa: ", material_params['kappa'])

    print('mesh parameters: time steps = %d, num_cells = %d, num_vertices = %d'
          % (test_params["time_steps"], mesh.num_cells(), mesh.num_vertices()))
    print('h_max = ', mesh.hmax())
    print('h_min = ', mesh.hmin())
    print("u_approx_order = ", test_params["u_approx_order"])
    print("p_approx_order = ", test_params["p_approx_order"])
    print("DOFs of u = ", len(func_spaces["V"].dofmap().dofs()))
    print("DOFs of p = ", len(func_spaces["Q"].dofmap().dofs()))
    print("DOFs of tau = ", len(func_spaces["H_Div"].dofmap().dofs()))
    print("DOFs of y = ", len(func_spaces["H_div"].dofmap().dofs()))

    print("iter_num = ", test_params["iter_num"])

def print_biot_iterative_coupling_parameters(beta, L, L_opt, q, q_opt):
    print('beta = ', float(beta))
    print('L = ', float(L))
    print('L_opt = ', float(L_opt))
    print('q = ', float(q))
    print('q_opt = ', float(q_opt))

# Outputting to the consol the results of the interation procedure
def output_errors_wrt_iterations(ep_array, epl2_array, eu_array, eudiv_array, ITER_NUM):

    print("")
    print("iter # & ||| e_p |||^2   & || e_p ||^2_{L2}  & ||| e_u |||^2   & || div e_u ||^2  ")
    print("----------------------------------------------------------------------------------")

    for i in range(1, ITER_NUM):
        print("    %2d & %10.4e & %10.4e & %10.4e  & %10.4e " % (i, ep_array[i], epl2_array[i], eu_array[i], eudiv_array[i]))
    print("")

# Outputting to the consol the results of the interation procedure
def output_errors_and_estimates_wrt_iterations(ep_array, epl2_array, eu_array, eudiv_array,
                                               maj_ep_array, maj_epl2_array, maj_eu_array, maj_eudiv_array,
                                               maj_pi_it, maj_ui_it,
                                               ITER_NUM, norm_p, norm_u):

    print("")
    print("iter # & ||| e_p |||^2   M^2_p     & || e_p ||^2_{L2}    M^2_{L2}  & ||| e_u |||^2   M^2_{u}    & || div e_u ||^2    M^2_{divu}")
    print("--------------------------------------------------------------------------------------------------------------------------------------")
    #for i in range(ITER_NUM-2, ITER_NUM):
    for i in range(0, ITER_NUM):
        print("    %2d & %10.4e &    %10.4e & %14.4e  &   %10.4e & %10.4e &    %10.4e  & %14.4e   &  %10.4e \\\\"
              % (i+1, ep_array[i] / norm_p, maj_ep_array[i] / norm_p,
                 epl2_array[i] / norm_p, maj_epl2_array[i] / norm_p,
                 eu_array[i] / norm_u, maj_eu_array[i] / norm_u,
                 eudiv_array[i] / norm_u, maj_eudiv_array[i] / norm_u))
    '''
    print("")
    print("final iter & ||| e_p |||^2   \tM^2_ph \t\tM^2_pi & ||| e_u |||^2    \tM^2_uh  \tM^2_ui ")
    print("--------------------------------------------------------------------------------------")
    for i in range(ITER_NUM-1, ITER_NUM):
        print("        %2d & %10.4e   %12.4e  %10.4e & %13.4e %12.4e %12.4e "
              % (i+1, ep_array[i], maj_ep_array[i], maj_pi_it[i],
                 eu_array[i], maj_eu_array[i],  maj_ui_it[i]))
    print("")
    '''

# -------------------------------------------------------------------------------------------------------------------- #
# Plotting functions
# -------------------------------------------------------------------------------------------------------------------- #
# Plotting the histogram
def plot_histogram(mesh, e_distr, m_distr, result_path):

    num_cells = mesh.num_cells()
    num_vert = mesh.num_vertices()
    cells = np.linspace(0, num_cells, num_cells)

    plt.figure()
    plt.plot(cells, e_distr, 'r.-')
    plt.plot(cells, m_distr, 'g.-')
    plt.xlabel('cells')
    plt.ylabel('e distr, m distr')
    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))
    plt.close('all')

# Plotting the scalar fucntion
def plot_function(f, mesh, dim, result_path):

    if dim == 1:
        plot_function_1d(mesh, f, result_path)
    elif dim == 2:
        #plot_function_2d(mesh, f, result_path)
        plot_function_3d(mesh, f, result_path)

def plot_function_1d(mesh, f, result_path):

    num_cells = mesh.num_cells()
    num_vert = mesh.num_vertices()

    fig = plt.gcf()
    plot(f)
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))

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

    plt.figure()
    plotting.plot(f, mode='color', cmap=cm.BrBG, edgecolors='k', mesh=mesh)

    # Create the triangulation
    #mesh_coordinates = mesh.coordinates().reshape((num_vert, d))
    #triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    #triangulation = tri.Triangulation(mesh_coordinates[:, 0],
    #                                  mesh_coordinates[:, 1],
    #                                  triangles)


    # Get the z values as face colors for each triangle(midpoint)
    '''
    plt.figure()
    zfaces = np.asarray([f(cell.midpoint()) for cell in dolfin.cells(mesh)])
    plt.tripcolor(tri.triangulation, facecolors=zfaces, edgecolors='k')
    fig = plt.gcf()

    fig.savefig(result_path + '-zfaces-cells-%d-vertices-%d.eps' %(num_cells, num_vert))
    '''
    # Get the z values for each vertex
    #plt.figure()
    #z = np.asarray([f(point) for point in mesh_coordinates])
    #plt.tripcolor(triangulation, z, cmap=cm.BrBG, edgecolors='k')
    #plt.colorbar()
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
    triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)
    z = np.asarray([f(point) for point in mesh_coordinates])

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


# Plotting the vector fucntion
def plot_vector_function(f, mesh, dim, result_path):

    if dim == 2:
        plot_vector_function_2d(mesh, f, result_path)

def plot_vector_function_2d(mesh, f, result_path):
    plt.figure()
    plotting.plot(f, cmap=cm.BrBG, mesh=mesh)
    # plt.colorbar(p)
    plt.savefig(result_path + "-mode.eps")

    plt.figure()
    plotting.plot(f, cmap=cm.BrBG, mode="displacement", mesh=mesh)
    # plt.colorbar()
    plt.savefig(result_path + "-mode-glyphs.eps")
    #plt.close('all')

# Plotting the mesh
def plot_mesh(mesh, result_path):

    plt.figure()
    #plt.hold(True)

    # change the font
    matplotlib.rcParams.update({'font.size': 10,
                                'font.family': 'serif'})

    n_vert = mesh.num_vertices()
    n_cells = mesh.num_cells()
    d = mesh.geometry().dim()

    plotting.plot(mesh)
    # Create the triangulation
    #mesh_coordinates = mesh.coordinates().reshape((n_vert, d))
    #triangles = np.asarray([cell.entities(0) for cell in cells(mesh)])
    #triangulation = tri.Triangulation(mesh_coordinates[:, 0],
    #                                  mesh_coordinates[:, 1],
    #                                  triangles)

    # Plot the mesh
    #plt.triplot(triangulation)
    plt.xlabel("x")
    plt.ylabel("y")

    # Saving the mesh into the file
    fig = plt.gcf()
    fig.savefig(result_path + '-cells-%d-vertices-%d' % (n_cells, n_vert) + ".eps")
    plt.close('all')


# Plotting the carpet with error and majorant distribution
def plot_carpet_2d(mesh, f, result_path):
    num_vert = mesh.num_vertices()
    num_cells = mesh.num_cells()

    d = mesh.geometry().dim()

    # Create the triangulation
    mesh_coordinates = mesh.coordinates().reshape((num_vert, d))
    triangles = np.asarray([cell.entities(0) for cell in dolfin.cpp.mesh.cells(mesh)])
    triangulation = tri.Triangulation(mesh_coordinates[:, 0],
                                      mesh_coordinates[:, 1],
                                      triangles)

    fig = plt.figure()

    # add labels
    plt.figure()
    plt.xlabel('x')
    plt.ylabel('y')

    # change the font
    matplotlib.rcParams.update({'font.size': 10,
                                'font.family': 'serif'})


    z = np.asarray([f(point) for point in mesh_coordinates])
    plt.tripcolor(triangulation, z, cmap=cm.coolwarm, edgecolors='k')
    fig = plt.gcf()
    fig.savefig(result_path + "-cells-%d-vertices-%d.eps" % (num_cells, num_vert))
    plt.close('all')

def mesh2triang(mesh):
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

def mplot_cellfunction(cellfn):
    C = cellfn.array()
    tri = mesh2triang(cellfn.mesh())
    return plt.tripcolor(tri, facecolors=C)

def mplot_function(f):
    mesh = f.function_space().mesh()
    if (mesh.geometry().dim() != 2):
        raise AttributeError('Mesh must be 2D')
    # DG0 cellwise function
    if f.vector().size() == mesh.num_cells():
        C = f.vector().array()
        return plt.tripcolor(mesh2triang(mesh), C)
    # Scalar function, interpolated to vertices
    elif f.value_rank() == 0:
        C = f.compute_vertex_values(mesh)
        return plt.tripcolor(mesh2triang(mesh), C, shading='gouraud')
    # Vector function, interpolated to vertices
    elif f.value_rank() == 1:
        w0 = f.compute_vertex_values(mesh)
        if (len(w0) != 2*mesh.num_vertices()):
            raise AttributeError('Vector field must be 2D')
        X = mesh.coordinates()[:, 0]
        Y = mesh.coordinates()[:, 1]
        U = w0[:mesh.num_vertices()]
        V = w0[mesh.num_vertices():]
        return plt.quiver(X,Y,U,V)

# Plot a generic dolfin object (if supported)

def plot(obj):
    plt.gca().set_aspect('equal')

    if isinstance(obj, Function):
        return mplot_function(obj)
    elif isinstance(obj, CellFunctionSizet):
        return mplot_cellfunction(obj)
    elif isinstance(obj, CellFunctionDouble):
        return mplot_cellfunction(obj)
    elif isinstance(obj, CellFunctionInt):
        return mplot_cellfunction(obj)
    elif isinstance(obj, Mesh):
        print("CHECK!")

    if (obj.geometry().dim() != 2):
        raise AttributeError('Mesh must be 2D')
        return plt.triplot(mesh2triang(obj), color='#808080')

    raise AttributeError('Failed to plot %s'%type(obj))

def create_line_mesh(vertices):
    "Given list of vertex coordinate tuples, build and return a mesh of intervals."

    # Get dimensions
    gdim = 1 if isinstance(vertices[0], float) else len(vertices[0])
    tdim = 1

    # Automatic choice of cellname for simplices
    cellname = "interval"

    # Indirect error checking and determination of tdim via ufl
    ufl_cell = Cell(cellname, gdim)
    assert tdim == ufl_cell.topological_dimension()

    # Create mesh to return
    mesh = Mesh()

    # Open mesh in editor
    me = MeshEditor()
    me.open(mesh, cellname, tdim, gdim)

    # Add vertices to mesh
    nv = len(vertices)
    me.init_vertices(nv)
    if gdim == 1:
        for i, v in enumerate(vertices):
            me.add_vertex(i, v)
    else:
        for i, v in enumerate(vertices):
            me.add_vertex(i, *v)

    # Add cells to mesh
    me.init_cells(nv-1)
    for i in range(nv-1):
        c = (i, i+1)
        me.add_cell(i, *c)

    me.close()

    return mesh

def line1d():
    n = 100
    us = [i/float(n-1) for i in range(n)]
    vertices = [u**3 for u in us]
    return create_line_mesh(vertices)

def line2d():
    n = 100
    us = [i/float(n-1) for i in range(n)]
    vertices = [(math.cos(1.5 * DOLFIN_PI * u),
                 math.sin(1.5 * DOLFIN_PI * u))
                 for u in us]
    return create_line_mesh(vertices)

def line3d():
    n = 100
    us = [i/float(n-1) for i in range(n)]
    vertices = [(math.cos(4.0 * DOLFIN_PI * u),
                 math.sin(4.0 * DOLFIN_PI * u),
                 2.0*u)
                 for u in us]
    return create_line_mesh(vertices)