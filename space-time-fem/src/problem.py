from dolfin.cpp.mesh import MeshFunction

__author__ = 'svetlana'

from dolfin import *
from numpy import info

expression_degree = 6
INFTY = 1e16

dirichlet_bc_marker = 1
neumann_bc_marker = 2
init_bc_marker = 3
do_nothing_bc_marker = 4

def SpaceDiv(f, dim):
    # Define differential operators for the 1D case
    if dim > 1:
        div_f = div(f)
    else:
        div_f = Dx(f, 0)
    return div_f

# define Div, Grad, Dt for space-time cylinder where time is the extra dimension
def Grad(f, dim):

    if dim == 2:
        # Define differential operator for the 1D case for space
        #space_matrix = as_matrix(((1, 0), (0, 0)))
        grad_f = f.dx(0)
    elif dim == 3:
        space_matrix = as_matrix(((1, 0, 0), (0, 1, 0), (0, 0, 0)))
        # Define differential operator for the 2D case for space
        grad_f = space_matrix * grad(f)
    #plotting.plot(grad(f))
    #grad_f = space_matrix * grad(f)

    return grad_f

def Laplas(f, dim):

    # Define Laplas for the different dimensions case for space
    if dim == 2:
        laplas_f = f.dx(0).dx(0)
    elif dim == 3:
        laplas_f = f.dx(0).dx(0) + f.dx(1).dx(1)

    return laplas_f


def NablaGrad(f, dim):
    if dim == 2:
        # Define differential operator for the 1D case for space
        # space_matrix = as_matrix(((1, 0), (0, 0)))
        grad_f = f.dx(0)
    elif dim == 3:
        space_matrix = as_matrix(((1, 0, 0), (0, 1, 0), (0, 0, 0)))
        # Define differential operator for the 2D case for space
        grad_f = nabla_grad(space_matrix * f)

    return grad_f

def D_t(f, dim):
    # Define time derivative for 'dim'-dimensional space
    return f.dx(dim-1)

def Div(f, dim):
    if dim == 2:
        div_f = f.dx(0)
    elif dim == 3:
        space_matrix = as_matrix(((1, 0, 0), (0, 1, 0), (0, 0, 0)))
        div_f = SpaceDiv(space_matrix * f, dim)
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

    if dim == 2:
        H_div = FunctionSpace(mesh, 'Lagrange', y_degree)
    elif dim == 3:
        H_div = VectorFunctionSpace(mesh, 'Lagrange', y_degree)

    return V, VV, V_exact, VV_exact, H_div

def output_problem_characteristics(test_num, u_expr, grad_u_expr, f_expr,
                                   domain, dim,
                                   num_cells, num_vertices, space_approx_deg):

    print "--------------------------------------------------------------------------------"
    print "test: ", test_num
    print "domain: ", domain
    print "u = ", u_expr
    print "f = ", f_expr

    print 'mesh parameters: num_cells = %d, num_vertices = %d' % (num_cells, num_vertices)
    print "space func approx_deg = ", space_approx_deg
    print "--------------------------------------------------------------------------------"


# Function to solve convection-reaction-diffusion
def solve_convection_reaction_diffusion(V, c_H, eps, f, A, min_eigs_A, lmbd, a,
                                        u0, uD, uN,
                                        mesh, boundary_facets,
                                        dim, test_num, test_params):
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

    # Define the stabilization parameter delta element-wise
    if test_params["STABILIZE"]:
        #fe = dolfin.TensorElement(family="Lagrange", cell=mesh.ufl_cell(), degree=3, quad_scheme="Default")
        fe = dolfin.FiniteElement("CG", cell=mesh.ufl_cell(), degree=3)
        #fe = dolfin.Element("Quadrature", cell=mesh.ufl_cell(), degree=3, quad_scheme="Default", dim=dim)
        delta = DeltaExpression(eps=eps, a=interpolate(Expression(str(c_H), degree=3, element=fe), V), mesh=mesh, element=fe)
        #delta = 0.5 * CellSize(mesh) / norm(c_H, 'linf')
    else:
        delta = 0.0

    # SUPG stabilization
    upwind = c_H * delta * D_t(v, dim) # with upwind
    #upwind = 0
    s_stab = delta * inner(c_H * D_t(u, dim)
                            - Div(A * NablaGrad(u, dim), dim)
                            + lmbd * u
                            + inner(a, NablaGrad(u, dim)),
                           upwind) * dx(domain=mesh)
    l_stab = delta * inner(f, upwind) * dx(domain=mesh)

    # Bilinear form
    a_stab = (c_H * inner(D_t(u, dim), v)
                + inner(A * NablaGrad(u, dim), NablaGrad(v, dim))
                + inner(lmbd * u, v)
                + inner(inner(a, NablaGrad(u, dim)), v)) * dx(domain=mesh) \
              + s_stab
    L_stab = (f * v) * dx(domain=mesh) \
              + (uN * v) * ds(neumann_bc_marker) \
              + l_stab

    '''
    theta = delta # depending on delta or just a constant
    test_func = v + c_H * theta * delta * D_t(v, dim)

    a_stab = (c_H * inner(D_t(u, dim), test_func),
              + inner(A * NablaGrad(u, dim), NablaGrad(test_func, dim))
              + inner(inner(a, NablaGrad(u, dim)), test_func)) * dx(domain=mesh)

    L_stab = (f * test_func) * dx(domain=mesh) \
             + (uN * test_func) * ds(neumann_bc_marker)
    '''
    # Define the unknown function
    u = Function(V)


    # Define boundary condition
    if test_num == 24 or test_num == 38 or test_num == 39 or test_num == 40 or test_num == 41 or test_num == 42: #  or test_num == 25
        bcs = [DirichletBC(V, u0, boundary_facets, init_bc_marker)]  # bottom of the cylinder with IC
    else:
        bcs = [DirichletBC(V, uD, boundary_facets, dirichlet_bc_marker),  # lateral surface Dirichlet BC
               DirichletBC(V, u0, boundary_facets, init_bc_marker)]  # bottom of the cylinder with IC
        # top of the cylinder is do-nothing condition

    # Solve the system generated by the variational equation a(u, v) = f(v)
    solve(a_stab == L_stab, u, bcs)

    return u, delta


def solve_parabolic_problem(V, c_H, eps, f, A, u0, uD, mesh, boundary_facets, dim, test_num):

    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)

    a = (c_H * inner(D_t(u, dim), v)
         + inner(A * NablaGrad(u, dim), NablaGrad(v, dim))) * dx(domain=mesh)
    L = f * v * dx(domain=mesh)

    u = Function(V)

    # Define boundary condition
    bcs = [DirichletBC(V, uD, boundary_facets, 0),  # lateral surface Dirichlet BC
           DirichletBC(V, u0, boundary_facets, 1)]  # bottom of the cylinder with IC
    """
    # Define boundary condition
    if test_num == 11 or test_num == 12 or test_num == 13:
        bcs = [DirichletBC(V, uD, boundary_facets, 0),
               DirichletBC(V, u0, boundary_facets, 1)]  # bottom with initial condition
    else:
        bcs = [DirichletBC(V, uD, boundary_facets, 0)]
    """
    solve(a == L, u, bcs)

    return u


def calculate_CF_of_domain(domain, dim):

    height = INFTY
    width  = INFTY
    length = INFTY

    if domain == "l-shape-domain" and dim == 2:

        height = 1.0
        width = 2.0

    elif domain == "1times1-minus-0times0" and dim == 2:

        height = 2.0
        width = 2.0

    elif domain == "unit-domain":

        if dim == 3:
            length = 1.0
            height = 1.0
            width = 1.0

        elif dim == 2:
            height = 1.0
            width = 1.0

        elif dim == 1:
            width = 1.0

    elif domain == "long-rectangle-domain":

        if dim == 2:
            width = 3.0

    elif domain == "rectangular-domain-1x2" and dim == 1:
        width = 2.0

    elif domain == "10byT-domain" and dim == 1:
        width = 2.0

    C_FD = 1.0 / DOLFIN_PI / sqrt(1.0 / height**2 + 1.0 / width**2 + 1.0 / length**2)
    print "C_FD = 1 / pi / (h_1^(-2) + ... + h_n^(-2)) = ", C_FD

    return C_FD

def calculate_trace_constant(test_num):

    # No Neumann part of the boundary
    C_Ntr = 0.0

    # Mark different parts of the domain for non-homogeneous Dirichlet BC
    if test_num == 48:
        # Caclulated in the matlab
        C_Ntr = 1.47

    if test_num == 38 and test_num == 39:
        # Caclulated in the matlab
        C_Ntr = 1.47
    return C_Ntr


# Class for delta-stabilization parameter (inherited from Expression)
class DeltaExpression(Expression):
    """
    :param eps: small parameter in fromt of diffusion operator, 0 < eps << 1
    :param a: convection field
    :param mesh: mesh
    :return: values of the delta function defined element-wise
    """

    def __init__(self, eps, a, mesh, element):
        #super(DeltaExpression, self).__init__(**kwargs)
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
        #a_norm = norm(self._a, 'linf')
        #a_norm = norm(self._a, 'L2')

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

# Interpolate scalar (dim = 1) and vector (dim > 1) functions
def interpolate_vector_function(a, dim, V_exact, VV_exact):
    """
    :param a: convection function
    :param dim: geometrical dimension
    :param V_exact: functional space of the exact function
    :param VV_exact: functional space of the exact vector function
    :return: interpolated vector fucntion
    """
    if dim == 2:
        #a_interpolated = a
        a_interpolated = interpolate(a, V_exact)
        #a_interpolated = project(a, V_exact)
    else:
        a_interpolated = interpolate(a, VV_exact)
        #a_interpolated = project(a, VV_exact)

    return a_interpolated
