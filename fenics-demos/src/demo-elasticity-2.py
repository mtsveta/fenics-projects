from dolfin import *

from dolfin.cpp.function import near
from dolfin.cpp.io import File, interactive
from dolfin.cpp.mesh import RectangleMesh, FacetFunction, Mesh

# Use UFLACS to speed-up assembly and limit quadrature degree
parameters['form_compiler']['representation'] = 'uflacs'
parameters['form_compiler']['optimize'] = True
parameters['form_compiler']['quadrature_degree'] = 4

def solve_elasticity(facet_function, E, nu, dt, T_end, output_dir):
    """Solves elasticity problem with Young modulus E, Poisson ration nu,
    timestep dt, until T_end and with output data going to output_dir.
    Geometry is defined by facet_function which also defines rest boundary
    by marker 1 and traction boundary by marker 2."""

    # Get mesh and prepare boundary measure
    mesh = facet_function.mesh()
    gdim = mesh.geometry().dim()
    ds = Measure("ds", subdomain_data=facet_function)

    # Build function space
    U = VectorFunctionSpace(mesh, "Lagrange", 1)
    P = FunctionSpace(mesh, "Lagrange", 1)
    W = MixedFunctionSpace([U, U, P])
    from dolfin.cpp.common import info
    info("Num DOFs %d" % W.dim())

    # Prepare BCs
    bcs = [DirichletBC(W.sub(i), gdim*(0.0,), facet_function, 1)
           for i in [0, 1]]

    # Define constitutive law
    def stress(u, p):
        """Returns 1st Piola-Kirchhoff stress and (local) mass balance
        for given u, p."""
        mu = Constant(E/(2.0*(1.0 + nu)))
        F = I + grad(u)
        J = det(F)
        B = F * F.T
        T = -p*I + mu*(B-I) # Cauchy stress
        S = J*T*inv(F).T # 1st Piola-Kirchhoff stress
        if nu == 0.5:
            # Incompressible
            pp = J-1.0
        else:
            # Compressible
            lmbd = Constant(E*nu/((1.0 + nu)*(1.0 - 2.0*nu)))
            pp = 1.0/lmbd*p + (J*J-1.0)
        return S, pp

    # Timestepping theta-method parameters
    q = Constant(0.5)
    dt = Constant(dt)

    # Unknowns, values at previous step and test functions
    w = Function(W)
    (u, v, p) = split(w)
    w0 = Function(W)
    (u0, v0, p0) = split(w0)
    (_u, _v, _p) = TestFunctions(W)

    I = Identity(W.mesh().geometry().dim())

    # Balance of momentum
    S, pp = stress(u, p)
    S0, pp0 = stress(u0, p0)
    F1 = (1.0/dt)*inner(u-u0, _u)*dx \
       - ( q*inner(v, _u)*dx + (1.0-q)*inner(v0, _u)*dx )
    F2a = inner(S, grad(_v))*dx + pp*_p*dx
    F2b = inner(S0, grad(_v))*dx + pp0*_p*dx
    F2 = (1.0/dt)*inner(v-v0, _v)*dx + q*F2a + (1.0-q)*F2b

    # Traction at boundary
    F = I + grad(u)
    bF_magnitude = Constant(0.0)
    bF_direction = {2: Constant((0.0, 1.0)), 3: Constant((0.0, 0.0, 1.0))}[gdim]
    bF = det(F)*dot(inv(F).T, bF_magnitude*bF_direction)
    FF = inner(bF, _v)*ds(2)

    # Whole system and its Jacobian
    F = F1 + F2 + FF
    J = derivative(F, w)

    # Initialize solver
    problem = NonlinearVariationalProblem(F, w, bcs=bcs, J=J)
    solver = NonlinearVariationalSolver(problem)
    solver.parameters['newton_solver']['relative_tolerance'] = 1e-6
    solver.parameters['newton_solver']['linear_solver'] = 'mumps'

    # Extract solution components
    (u, v, p) = w.split()
    u.rename("u", "displacement")
    v.rename("v", "velocity")
    p.rename("p", "pressure")

    # Create files for storing solution
    vfile = File("%s/velo.xdmf" % output_dir)
    ufile = File("%s/disp.xdmf" % output_dir)
    pfile = File("%s/pres.xdmf" % output_dir)

    # Prepare plot window
    plt = plot(u, mode="displacement", interactive=False, wireframe=True)

    # Time-stepping loop
    t = 0.0
    while t <= T_end:
        print "Time: %g" % t
        t += float(dt)

        # Increase traction
        bF_magnitude.assign(100.0*t)

        # Prepare to solve and solve
        w0.assign(w)
        solver.solve()

        # Store solution to files and plot
        ufile << (u, t)
        vfile << (v, t)
        pfile << (p, t)
        plt.plot(u)


def geometry_2d(length):
    """Prepares 2D geometry. Returns facet function with 1, 2 on parts of
    the boundary."""
    n = 4
    x0 = 0.0
    x1 = x0 + length
    y0 = 0.0
    y1 = 1.0
    mesh = RectangleMesh(x0, y0, x1, y1, int((x1-x0)*n), int((y1-y0)*n), 'crossed')
    boundary_parts = FacetFunction('size_t', mesh)
    left  = AutoSubDomain(lambda x: near(x[0], x0))
    right = AutoSubDomain(lambda x: near(x[0], x1))
    left .mark(boundary_parts, 1)
    right.mark(boundary_parts, 2)
    boundary_parts._mesh = mesh # Workaround issue #467
    return boundary_parts


def geometry_3d():
    """Prepares 3D geometry. Returns facet function with 1, 2 on parts of
    the boundary."""
    mesh = Mesh('lego_beam.xml')
    gdim = mesh.geometry().dim()
    x0 = mesh.coordinates()[:, 0].min()
    x1 = mesh.coordinates()[:, 0].max()
    boundary_parts = FacetFunction('size_t', mesh)
    left  = AutoSubDomain(lambda x: near(x[0], x0))
    right = AutoSubDomain(lambda x: near(x[0], x1))
    left .mark(boundary_parts, 1)
    right.mark(boundary_parts, 2)
    boundary_parts._mesh = mesh # Workaround issue #467
    return boundary_parts


if __name__ == '__main__':

    solve_elasticity(geometry_2d(20.0), 1e5, 0.3, 0.25, 5.0, 'results_2d_comp')
    solve_elasticity(geometry_2d(20.0), 1e5, 0.5, 0.25, 5.0, 'results_2d_incomp')
    solve_elasticity(geometry_2d(80.0), 1e5, 0.3, 0.25, 5.0, 'results_2d_long_comp')
    solve_elasticity(geometry_3d(),     1e5, 0.3, 0.50, 5.0, 'results_3d_comp')
    interactive()