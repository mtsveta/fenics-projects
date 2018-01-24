__author__ = 'svetlana'

from dolfin import *
from dolfin.cpp.mesh import MeshFunction, cells, CellFunction
from ufl.tensoralgebra import Inverse
from numpy.linalg import inv
import numpy, postprocess, math

expression_degree = 4

# Generalized inverse operator
def inverse(A, A_expr, dim):
    """
    :param A: 1d-scalar or nd-matrix to inverse
    :param dim: dimension of the problem
    :return: inverted A depending on dim
    """

    if dim == 1:
        invA = 1 / A
    else:
        A_mtx = numpy.matrix(A_expr)
        inv_A_mtx = inv(A_mtx)
        invA = as_matrix(inv_A_mtx.tolist())

    return invA

# Generalized gradient operator
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

# Generalized divergence operator
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

# Function to generate functional spaces
def functional_spaces(mesh, test_params, dim):
    """
    :param mesh: mesh
    :param v_degree: degree of basis ussed to approximate approximation v
    :param y_degree: degree of basis ussed to approximate flux y
    :param dim: dimension of the problem
    :return: constructed functional spaces
             V approximation space
             VV approximation vector-space
             V_exact high order approximation space
             VV_exact high order approximation vector-space
             H_div space for fluxes
    """
    v_degree = test_params['v_approx_order']
    y_degree = test_params['flux_approx_order']
    # Order of the space of the exact solutions
    v_exact_degree = test_params["v_exact_approx_order"]

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

# Output description of the problem
def output_problem_characteristics(test_num, u_expr, grad_u_expr, f_expr, A_expr, lambda_expr, a_expr, uD_expr,
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
    print "a = ", a_expr
    print "uD = ", uD_expr
    print "u = ", u_expr
    print "grad_u = ", grad_u_expr


    print 'mesh parameters: num_cells = %d, num_vertices = %d' % (num_cells, num_vertices)
    print "space func approx_deg = ", space_approx_deg
    print "flux approx_deg = ", flux_approx_deg
    print "dim", dim

# Calculate Friedrichs constant of thegdiven domian
def calculate_CF_of_domain(domain, dim):
    """
    :param domain: domain name
    :param dim: dimension of the problem
    :return: C_FD Friedrichs constant
    """
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

    print "C_FD = 1 / pi / (h_1^(-2) + ... + h_n^(-2)) = ", C_FD

    return C_FD

# Class for delta-stabilization parameter (inherited from Expression)
class DeltaExpression(Expression):
    """
    :param eps: small parameter in fromt of diffusion operator, 0 < eps << 1
    :param a: convection field
    :param mesh: mesh
    :return: values of the delta function defined element-wise
    """

    def __init__(self, eps, a, mesh, element):
        #super(DeltaExpression, self).__init__(degree=expression_degree)
        self._eps = eps
        self._a = a
        self._mesh = mesh
        self._element = element

    # Redefine eval_cell function
    def eval_cell(self, value, x, ufc_cell):
        """
        :param value: values of the delta-function element-wise
        :param x:
        :param ufc_cell: cell of the mesh
        :return: values of the delta function defined element-wise
        """

        # Obtain the cell based on its index in the mesh
        cell = Cell(self._mesh, ufc_cell.index)
        # Get the diameter of the cell
        h = cell.h()
        # Calculate the || a ||_inf norm
        a_norm = norm(self._a.vector(), 'linf')

        # Calculate Peclet number
        Pe = h * a_norm / 2 / self._eps

        # Defined auxilary constants
        C_0 = 0.1           # 0.1 in L.Chen paper
        C_1 = 0.0    # 0.0 in L.Chen paper

        # Defined stabilized parameter depending on Peclet number
        if Pe > 1:
            value[0] = C_0 * h / a_norm # convection dominated
        elif Pe <= 1:
            value[0] = C_1 * (h)**2 / self._eps # diffusion dominated

# Class to detecting the smoothness (inherited from Expression)
class SmoothnessCriterion(Expression):
    """
    :param v: function
    :param mesh: mesh
    :return: values of the delta function defined element-wise
    """

    def __init__(self, v, mesh, expression_degree):
        #super(SmoothnessCriterion, self).__init__(degree=expression_degree)
        #super.__init__(degree=expression_degree)
        self._v = v
        self._mesh = mesh

    # Redefine eval_cell function
    def eval_cell(self, value, x, ufc_cell):
        """
        :param value: values of the delta-function element-wise
        :param x:
        :param ufc_cell: cell of the mesh
        :return: values of the delta function defined element-wise
        """

        # Obtain the cell based on its index in the mesh
        cell = Cell(self._mesh, ufc_cell.index)
        # Get the diameter of the cell
        h = cell.h()
        dim = cell.num_vertices() - 1

        # Calculate the required norms
        v = self._v
        v_inf_norm    = l_inf_norm(v, cell.get_vertex_coordinates().reshape(dim + 1, dim), dim)
        v_h1_seminorm = assemble_local(inner(grad(v), grad(v)) * dx, cell)
        v_l2_norm     = assemble_local(v * v * dx, cell)

        #print "v_inf_norm**2", v_inf_norm**2
        #print "v_h1_seminorm", v_h1_seminorm
        #print "v_l2_norm", v_l2_norm
        # coth = 1 / tanh
        value[0] = v_inf_norm**2 / (1/math.tanh(1) * (h * v_h1_seminorm + v_l2_norm / h))

def coth(x):
    return 1./math.tanh(x)

# Function to solve convection-reaction-diffusion
def solve_convection_reaction_diffusion(V, f, A, lambda_l, lmbd, a,
                                        boundary, uD, uN,
                                        mesh, boundary_facets,
                                        dim, test_num,
                                        test_params):
    """
    :param V: functional space of the approximation
    :param f: right-hand side
    :param A: diffusion matrix
    :param lambda_l: minimal eigenvalue of matrix A
    :param lmbd: reaction function
    :param a: convection function
    :param uD: Dirichlet boundary function
    :param uD_boundary: definition of the Dirichlet boundary
    :param mesh: mesh-discretization of the domain
    :param dim: geometrical dimension of the domain
    :return: solution u and stabilization element-wise fucntion delta
    """

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)

    #U = u.vector()

    # Define the stabilization parameter delta element-wise
    if test_params["STABILIZE"]:

        fe = dolfin.FiniteElement("CG", cell=mesh.ufl_cell(), degree=3)
        # fe = dolfin.Element("Quadrature", cell=mesh.ufl_cell(), degree=3, quad_scheme="Default", dim=dim)
        delta = DeltaExpression(eps=lambda_l,
                                a=a,
                                mesh=mesh,
                                element=fe)
        # delta = 0.5 * CellSize(mesh) / norm(c_H, 'linf')
        #delta = DeltaExpression(eps=lambda_l, a=a, mesh=mesh, expression_degree=test_params["expression_degree"])
        delta = 0.5 * CellSize(mesh) / norm(a.vector(), 'linf')
    else:
        delta = 0.0

    # If v \in P1: divAgradu = 0
    if test_params["v_approx_order"] == 1:
        divAgradu = 0.0
    else:
        # Define bilinear form a(u, v) and linear functional f(v)
        if dim == 1:
            divAgradu = - A * Div(Grad(u, dim), dim)
        else:
            divAgradu = - Div(A * Grad(u, dim), dim)

    # SUPG stabilization
    #"""
    s_stab = delta * (inner(divAgradu + lmbd * u + inner(a, Grad(u, dim)),
                      inner(a, Grad(v, dim)))) * dx(domain=mesh)
    l_stab = delta * (inner(f, inner(a, Grad(v, dim)))) * dx(domain=mesh)
    #"""
    # GLS stabilization
    """
    s_stab = delta * (inner(divAgradu + lmbd * u + inner(a, Grad(u, dim)),
                      inner(a, Grad(v, dim)) + lmbd * v)) * dx(domain=mesh)
    l_stab = delta * (inner(f, inner(a, Grad(v, dim)) + lmbd * v)) * dx(domain=mesh)
    """

    # CIP stabilization
    """
    n = FacetNormal(mesh)
    # Check this jump fucntion
    s_stab = delta * (inner(inner(a, (Grad(jump(v, n), dim))), inner(a, Grad(v, dim)))) * dS(domain=mesh),
    l_stab = 0
    """

    # CIP stabilization
    """
    n = Normal(mesh)
    s_stab = delta * (inner(inner(n, Grad(v, dim)),
                      ) * dS(domain=mesh),
    l_stab = 0
    """

    # Bilinear form
    a_stab = (inner(A * Grad(u, dim), Grad(v, dim)) + lmbd * inner(u, v) + inner(inner(a, Grad(u, dim)), v)) * dx(domain=mesh) \
              + s_stab

    if test_num == 48:
        L_stab = (f * v) * dx(domain=mesh) + (uN * v) * ds(1) \
                  + l_stab
    else:
        L_stab = (f * v) * dx(domain=mesh) \
                  + l_stab
    # Define the unknown function
    u = Function(V)

    # Define boundary condition
    if test_num == 45:
        # Dirichlet boudary is devided intto south (bottom) part and rest of it
        bc = [DirichletBC(V, uD[0], boundary_facets, 0),
              DirichletBC(V, uD[1], boundary_facets, 1)]
    elif test_num == 48:
        # Dirichlet boundary is all parts of boudary except the east(right) part of it
        bc = DirichletBC(V, uD, boundary_facets, 0)
    else:
        bc = DirichletBC(V, uD, boundary)

    problem = LinearVariationalProblem(a_stab, L_stab, u, bc)
    solver = LinearVariationalSolver(problem)
    solver.parameters.linear_solver = 'gmres'
    solver.parameters.preconditioner = 'ilu'
    prm = solver.parameters.krylov_solver
    prm.absolute_tolerance = 1e-10
    prm.relative_tolerance = 1e-6
    prm.maximum_iterations = 1000

    # Solve the system generated by the variational equation a(u, v) = f(v)
    if len(u.vector()) <= 1e4:
        solve(a_stab == L_stab, u, bcs=bc)
    else:
        solver.solve()

        #solve(a_stab == L_stab, u, bcs=bc,
        #  solver_parameters={"linear_solver": "gmres", "preconditioner": "ilu"},
        #  form_compiler_parameters={"optimize": True})

    """
    for c in cells(mesh):
        # vert: num of vertices
        # horiz: dimension
        verts = c.get_vertex_coordinates().reshape(dim + 1, dim)
        u_inf_norm = l_inf_norm(u, verts, dim)

        h = c.h()
        #print "h ", h
        denominator = h / 2 * assemble_local(inner(grad(u), grad(u)) * dx, c) + 2 / h * assemble_local(u * u * dx, c)
        #print "denominator", denominator
        F_c = u_inf_norm**2 / denominator
        print "F_c: ", F_c
        num = c.index()
        #smooth_criterion[c] = F_c
    """
    return u, delta

def l_inf_norm(v, coords, dim):
    abs_vals = postprocess.allocate_array(dim + 1)
    for i in range(0, dim+1):
        abs_vals[i] = abs(v(coords[i, :]))
    return numpy.max(abs_vals)

def smoothness_distribution(u, mesh, dim):
    cell_num = mesh.num_cells()
    smoothness = postprocess.allocate_array(cell_num)

    for c in cells(mesh):
        # vert: num of vertices
        # horiz: dimension
        verts = c.get_vertex_coordinates().reshape(dim + 1, dim)
        u_inf_norm = l_inf_norm(u, verts, dim)

        h = c.h()
        denominator = h / 2 * assemble_local(inner(grad(u), grad(u)) * dx, c) + 2 / h * assemble_local(u * u * dx, c)
        F_c = u_inf_norm ** 2 / denominator
        smoothness[c.index()] = F_c
    return smoothness

# Interpolate scalar (dim = 1) and vector (dim > 1) functions
def interpolate_vector_function(a, dim, V_exact, VV_exact):
    """
    :param a: convection function
    :param dim: geometrical dimension
    :param V_exact: functional space of the exact function
    :param VV_exact: functional space of the exact vector function
    :return: interpolated vector fucntion
    """
    if dim == 1:
        a_interpolated = interpolate(a, V_exact)
    else:
        a_interpolated = interpolate(a, VV_exact)

    return a_interpolated

# Function that construct diffusion matrix as the cell-wise function
def construct_from_mesh_functions(dim, A_expr, mesh):
    """
    :param dim: geometrical dimension
    :param A_expr: expression defining diffusion matrix
    :param mesh: mesh-discretization of the domain
    :return: cell-wise function
    """

    # Define scalar cell-wise function for dim = 1
    a00 = MeshFunction("double", mesh, dim)
    # Define function expression from the problem data (depending on how many of the conditions are discussed )
    A1 = A_expr[0]
    A2 = A_expr[1]

    # Case for dimension higher then one (symmetric case)
    if dim >= 2:
        a01 = MeshFunction("double", mesh, dim)
        a11 = MeshFunction("double", mesh, dim)
    # Case for dimension higher then two (symmetric case)
    if dim >= 3:
        a02 = MeshFunction("double", mesh, dim)
        a12 = MeshFunction("double", mesh, dim)
        a22 = MeshFunction("double", mesh, dim)

    for cell in cells(mesh):
        if cell.midpoint().x() < 0.5: # this condition checks whethe the elements on the left part of the domain
            A = A_expr[0]
        else:
            A = A_expr[1]

        a00[cell] = A[0][0]
        if dim >= 2:
            a11[cell] = A[1][1]
            a01[cell] = A[0][1]
        if dim >= 3:
            a02[cell] = A[0][2]
            a12[cell] = A[1][2]
            a22[cell] = A[2][2]

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
    # Define expression via C++ code
    a = Expression(cppcode=conductivity_code)
    a.a00 = a00
    if dim >= 2:
        a.a01 = a01
        a.a11 = a11
    if dim >= 3:
        a.a02 = a02
        a.a12 = a12
        a.a22 = a22
    # Define matrix depending on the dimension
    if dim == 1:
        A = a[0]
    elif dim == 2:
        # A = |a[0] a[1]|
        #     |a[1] a[2]|
        A = as_matrix(((a[0], a[1]), (a[1], a[2])))
    elif dim == 3:
        # A = |a[0] a[1] a[2]|
        #     |a[1] a[3] a[4]|
        #     |a[2] a[4] a[5]|
        A = as_matrix(((a[0], a[1], a[2]), (a[1], a[3], a[4]), (a[2], a[3], a[5])))

    return A

# Contruction of system matrixes
def construct_stiffness_and_mass_matrices(u, v, lmbd):
    """
    :param u: function
    :param v: test function
    :param v: test function
    :param lmbda: reaction function
    :return: stiffness and mass matrixes
    """

    # Variational forms of stiffness (K), mass (M) matrices
    a_K = inner(nabla_grad(u), nabla_grad(v)) * dx
    a_M = u * v * dx
    a_R = lmbd * u * v * dx

    K = assemble(a_K)
    M = assemble(a_M)
    R = assemble(a_R)

    return K, M, R

