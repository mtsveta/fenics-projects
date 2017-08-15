from dolfin.cpp.function import near

__author__ = 'svetlana'

from dolfin import *
from ufl.tensoralgebra import Inverse
from dolfin.cpp.mesh import MeshFunction, cells, RectangleMesh, Point, FacetFunction
import mshr

# Identity tensor
def I_tensor(mesh):
    return Identity(mesh.geometry().dim())

# Strain tensor
def strain_tensor(v, dim):
    return sym(Grad(v, dim))

#Construct 2D geometry
def geometry_2d(domain_params, test_num, res):

    if domain_params["domain_type"] == "rectangular-domain":

        # Define coordinates of a rectangular domain
        x0 = 0.0
        x1 = x0 + domain_params["l_x"]
        y0 = 0.0
        y1 = y0 + domain_params["l_y"]

        # irregular mesh
        mesh_type = 'crossed'
        mesh = RectangleMesh(Point(x0, y0), Point(x1, y1), 20, 5, mesh_type)

        # regular mesh
        #resolution = res
        #mesh = mshr.generate_mesh(mshr.Rectangle(Point(x0, y0), Point(x1, y1)), resolution)

        # Get boundary facets
        boundary_parts = FacetFunction('size_t', mesh)
        boundary_parts.set_all(0)

        if test_num == 1:

            # Creat parts of the boundary as subdomains
            left  = AutoSubDomain(lambda x: near(x[0], x0))
            right = AutoSubDomain(lambda x: near(x[0], x1))
            bottom  = AutoSubDomain(lambda x: near(x[1], y0))
            top = AutoSubDomain(lambda x: near(x[1], y1))

            # Mark part of the domain with Neumann condition
            # dy default all Dirichlet part is marked as 0
            left.mark(boundary_parts, 1)
            right.mark(boundary_parts, 2)
            bottom .mark(boundary_parts, 3)
            top.mark(boundary_parts, 4)

        elif test_num == 3:

            # Creat parts of the boundary as subdomains
            gamma_1 = AutoSubDomain(lambda x: near(x[1], y1) and x[0] >= 0.25 and x[0] <= 0.75)
            gamma_1.mark(boundary_parts, 1)

        #boundary_parts._mesh = mesh
        #plot(mesh, interactive=True)

    elif domain_params["domain_type"] == "rectangular-with-obstacle":

        # Define basic rectangulars
        big_rect = mshr.Rectangle(Point(-2.0, 0.0), Point(2.0, 4.0))
        long_rect = mshr.Rectangle(Point(0.0, 0.0), Point(1.0, 2.0))

        # Define result geometry
        geometry = big_rect - long_rect

        # Build mesh
        resolution = 50
        mesh = mshr.generate_mesh(geometry, resolution)

    # we can send the facet-funtion, which contain the mesh
    return boundary_parts

def construct_stiffness_and_mass_matrices(u, v, A, dim):
    # Components of stiffness (K) & mass (M) matrices
    a_K = inner(A * Grad(u, dim), Grad(v, dim)) * dx
    a_M = u * v * dx

    K = assemble(a_K)
    M = assemble(a_M)
    return K, M

def output_problem_characteristics(test_num, f_expr, A_expr, lambda_expr, uD_expr,
                                   domain, dim, num_cells, num_vertices, space_approx_deg, flux_approx_deg):

    print " "
    print "%-------------------------"
    print "% Problem characteristics"
    print "%-------------------------"

    print "test: ", test_num
    print "domain: ", domain
    print "f = ", f_expr
    print "lambda = ", lambda_expr
    print "A = ", A_expr
    print "uD = ", uD_expr

    print 'mesh parameters: num_cells = %d, num_vertices = %d' % (num_cells, num_vertices)
    print "space func approx_deg = ", space_approx_deg
    print "flux approx_deg = ", flux_approx_deg


def calculate_CF_of_domain(domain_params):

    if domain_params["domain_type"] == "rectangular-domain" and domain_params["gdim"] == 2:

        height = domain_params["l_y"]
        width = domain_params["l_x"]
        C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / width**2 + 1.0 / height**2)

    print "Friedrichs' constant C_FD = 1 / pi / (h_1^(-2) + ... + h_n^(-2)) = ", C_FD

    return C_FD

def Grad(f, dim):
    """
    :param f: function
    :param dim: dimension of the problem
    :return: gradient of the function depending on dim
    """
    if dim > 1:
        grad_f = grad(f)
    else:
        # Define grad operator for the 1D case
       grad_f = Dx(f, 0)
    return grad_f

def Div(f, dim):
    """
    :param f: function
    :param dim: dimension of the problem
    :return: divergence of the function depending on dim
    """
    if dim > 1:
        div_f = div(f)
    else:
        # Define divergent operator for the 1D case
        div_f = Dx(f, 0)
    return div_f

def inverse(A, dim):
    if dim == 1:
        invA = 1 / A
    else:
        invA = Inverse(A)
    return invA

def get_next_step_solution(v_k, v_k1, K, tau, M, f_k, f_k1, bc, time_discretization_tag):

    if time_discretization_tag == 'implicit':

        # IMPLICIT time-discretization scheme
        # Compose the matrix of system of linear equations
        A = M + tau * K
        # Update the right-hand side of the system of linear equation
        b = M * v_k.vector() + tau * M * f_k1.vector()

    elif time_discretization_tag == 'explicit':
        # EXPLICIT time-discretization scheme
        # Update right hand side f_{k + 1}
        # Compose the matrix of system of linear equations
        A = M
        # Update the right-hand side of the system of linear equation
        b = M * v_k.vector() + tau * K * v_k.vector() + tau * M * f_k.vector()


    # Solve system of linear equations
    '''
    print "M = ", M.array()
    print "K = ", K.array()
    print "A = ", A.array()
    print "b = ", b.array()
    list_linear_solver_methods()
    list_krylov_solver_preconditioners()
    '''
    bc.apply(A, b)
    solve(A, v_k1.vector(), b, "petsc")

    return v_k1


def construct_A_from_mesh_functions(dim, A_expr, mesh):

    if dim == 2:
        a01 = MeshFunction("double", mesh, dim)
        a11 = MeshFunction("double", mesh, dim)

    A1 = A_expr[0]
    A2 = A_expr[1]

    a00 = MeshFunction("double", mesh, dim)

    for cell in cells(mesh):
        if cell.midpoint().x() < 0.5:
            if dim == 2:
                a00[cell] = A1[0][0]
                a11[cell] = A1[1][1]
                a01[cell] = A1[0][1]
        else:
            if dim == 2:
                a00[cell] = A2[0][0]
                a11[cell] = A2[1][1]
                a01[cell] = A2[0][1]

    # Code for C++ evaluation of conductivity
    conductivity_code = """

    class Conductivity : public Expression
    {
    public:

      // Create expression with 3 components
      Conductivity() : Expression(3) {}

      // Function for evaluating expression on each cell
      void eval(Array<double>& values, const Array<double>& x, const ufc::cell& cell) const
      {
        const uint D = cell.topological_dimension;
        const uint cell_index = cell.index;
        values[0] = (*a00)[cell_index];
        values[1] = (*a01)[cell_index];
        values[2] = (*a11)[cell_index];
      }

      // The data stored in mesh functions
      std::shared_ptr<MeshFunction<double> > a00;
      std::shared_ptr<MeshFunction<double> > a01;
      std::shared_ptr<MeshFunction<double> > a11;

    };
    """

    a = Expression(cppcode=conductivity_code)
    a.a00 = a00
    a.a01 = a01
    a.a11 = a11
    A = as_matrix(((a[0], a[1]), (a[1], a[2])))

    return A
