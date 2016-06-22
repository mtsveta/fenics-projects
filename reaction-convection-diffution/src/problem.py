from ufl.tensors import as_scalar

__author__ = 'svetlana'

from dolfin import *
from numpy import info
from dolfin.cpp.mesh import MeshFunction, cells, CellFunction
from ufl.tensoralgebra import Inverse
import ufl


def inverse(A, dim):
    if dim == 1:
        invA = 1 / A
    else:
        invA = Inverse(A)
    return invA

def Grad(f, dim):
    """
    :param f: function
    :param dim: dimension of the problem
    :return: gradient of the function depending on dim
    """
    if dim > 1:
        div_f = grad(f)
    else:
        # Define grad operator for the 1D case
        div_f = Dx(f, 0)
    return div_f

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

# Contruction of system matrixes
def construct_stiffness_and_mass_matrices(u, v, lmbd):
    # Variational forms of stiffness (K) & mass (M) matrices
    a_K = inner(nabla_grad(u), nabla_grad(v)) * dx
    a_M = u * v * dx

    a_R = lmbd * u * v * dx

    K = assemble(a_K)
    M = assemble(a_M)
    R = assemble(a_R)

    return K, M, R


def functional_spaces(mesh, v_degree, y_degree, dim):

    # Order of the space of the exact solutions
    v_exact_degree = v_degree + 3

    # Define functional spaces
    V = FunctionSpace(mesh, 'Lagrange', v_degree)
    VV = VectorFunctionSpace(mesh, 'Lagrange', v_degree)

    V_exact = FunctionSpace(mesh, 'Lagrange', v_exact_degree)
    VV_exact = VectorFunctionSpace(mesh, 'Lagrange', v_exact_degree)

    if dim == 1:
        H_div = FunctionSpace(mesh, 'Lagrange', y_degree)
    else:
        H_div = FunctionSpace(mesh, 'RT', y_degree)


    return V, VV, V_exact, VV_exact, H_div

def output_problem_characteristics(test_num, u_expr, grad_u_expr, f_expr, A_expr, lambda_expr, a_expr, uD_expr,
                                   domain, dim, num_cells, num_vertices, space_approx_deg, flux_approx_deg):

    print " "
    print "%-------------------------"
    print "% Problem characteristics"
    print "%-------------------------"

    print "test: ", test_num
    print "domain: ", domain
    print "u = ", u_expr
    print "f = ", f_expr
    print "lambda = ", lambda_expr
    print "A = ", A_expr
    print "a = ", a_expr
    print "uD = ", uD_expr

    print 'mesh parameters: num_cells = %d, num_vertices = %d' % (num_cells, num_vertices)
    print "space func approx_deg = ", space_approx_deg
    print "flux approx_deg = ", flux_approx_deg


def calculate_CF_of_domain(domain, dim):

    if domain == "l-shape-domain" and dim == 3:

        height = 1
        width = 2
        length = 2
        C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height**2 + 1.0 / width**2 + 1.0 / length**2)

    elif domain == "1times1-minus-0times0" and dim == 3:

        height = 2
        width = 2
        length = 2

        C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height**2 + 1.0 / width**2 + 1.0 / length**2)

    elif domain == "unit-domain":

        height = 1
        width = 1
        length = 1

        if dim == 3:
            C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height**2 + 1.0 / width**2 + 1.0 / length**2)
        elif dim == 2:
            C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / width**2 + 1.0 / length**2)
        elif dim == 1:
            C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / length**2)

    elif domain == "l-shape-domain" and dim == 2:

        width = 2
        length = 2

        C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / width**2 + 1.0 / length**2)

    elif domain == "circle-domain" and dim == 2:

        width = 4
        length = 4

        C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / width**2 + 1.0 / length**2)

    print "domain: ", domain
    print "C_FD = 1 / pi / (h_1^(-2) + ... + h_n^(-2)) = ", C_FD

    return C_FD

def solve_convection_reaction_diffusion(V, f, A, lambda_l, lmbd, a, u0, u0_boundary, mesh, dim):

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)

    def delta_T(mesh, a, lambda_l):

        h = CellSize(mesh)
        C_1 = 0.2
        delta = C_1 * h * h / lambda_l
        return delta

    class DeltaExpression(Expression):
        def __init__(self, lambda_min, a, mesh):
            self._lamda_min = lambda_min
            self._a = a
            self._mesh = mesh


        def eval_cell(self, value, x, ufc_cell):
            cell = Cell(self._mesh, ufc_cell.index)
            h = cell.diameter()
            a_norm = norm(a.vector(), 'linf', mesh=mesh)
            Pe = h * a_norm / 2 / self._lamda_min
            C_0 = 1e-2
            C_1 = 1e-2
            if Pe > 1:
                value[0] = C_0 * h
            elif Pe <= 1:
                value[0] = C_1 * (h)**2 / self._lamda_min

    #delta = DeltaExpression(lambda_min=lambda_l, a=a, mesh=mesh)
    delta = 0

    a_pert = (inner(A * Grad(u, dim), Grad(v, dim)) + lmbd * inner(u, v) + inner(inner(a, Grad(u, dim)), v) +
              delta * (inner(Div(Grad(u, dim), dim) + lmbd * u + inner(a, Grad(u, dim)), inner(a, Grad(v, dim))))) * dx(domain=mesh)
    L_pert = (f * v +
              delta * (inner(f, inner(a, Grad(v, dim))))) * dx(domain=mesh)

    u = Function(V)

    # Define boundary condition
    bc = DirichletBC(V, u0, u0_boundary)

    #solve(a == L, u, bc)

    solve(a_pert == L_pert, u, bc)
    #plot(u)
    #print "delta", delta

    return u, delta

# Interpolate scalar and vector functions
def interpolate_vector_function(a, dim, V_exact, VV_exact):

    if dim == 1:
        a_interpolated = interpolate(a, V_exact)
    else:
        a_interpolated = interpolate(a, VV_exact)

    return a_interpolated

def construct_from_mesh_functions(dim, A_expr, mesh):

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

