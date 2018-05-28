# -*- coding: utf-8 -*-
# Copyright (C) 2008-2012 Joachim B. Haga and Fredrik Valdmanis
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Martin Sandve Alnæs 2008-2015
# Modified by Anders Logg 2008-2015
#
# First added:  2008-03-05
# Last changed: 2015-11-26

from __future__ import print_function
import os
from distutils.version import StrictVersion

import dolfin
import dolfin.cpp as cpp
import ufl
import numpy as np

__all__ = ['plot']

_meshfunction_types = (cpp.MeshFunction, cpp.MeshFunctionBool, cpp.MeshFunctionInt,
                      cpp.MeshFunctionDouble, cpp.MeshFunctionSizet)
_plottable_types = (cpp.Function,
                    cpp.Expression,
                    cpp.Mesh,
                    cpp.MultiMesh,
                    cpp.DirichletBC) + _meshfunction_types

def _has_matplotlib():
    try:
        import matplotlib
    except ImportError:
        return False

    # Switch to Agg backend if DISPLAY not set
    if not os.environ.get("DISPLAY"):
        cpp.log(cpp.PROGRESS, "Environment variable DISPLAY not set. Switching "
                              "to 'Agg' matplotlib backend.")
        try:
            matplotlib.use("Agg")
        except ValueError as e:
            cpp.warning("Switching to 'Agg' backend failed with the message:")
            print('"%s"' % str(e))
            cpp.warning("Trying to continue...")

    return True

# Add plotting_backend to parameter system (don't need this in C++)
'''
_available_backends = []
if cpp.has_vtk():
    _available_backends.append("vtk")
if _has_matplotlib():
    _available_backends.append("matplotlib")
_available_backends.append("none")
cpp.parameters.add("plotting_backend", _available_backends[0],
        _available_backends)
'''
def mesh2triang(mesh):
    import matplotlib.tri as tri
    xy = mesh.coordinates()
    return tri.Triangulation(xy[:, 0], xy[:, 1], mesh.cells())

def mplot_mesh(ax, mesh, **kwargs):
    tdim = mesh.topology().dim()
    gdim = mesh.geometry().dim()
    if gdim == 2 and tdim == 2:
        color = kwargs.pop("color", '#808080')
        return ax.triplot(mesh2triang(mesh), color=color, **kwargs)
    elif gdim == 3 and tdim == 3:
        bmesh = dolfin.BoundaryMesh(mesh, "exterior", order=False)
        mplot_mesh(ax, bmesh, **kwargs)
    elif gdim == 3 and tdim == 2:
        xy = mesh.coordinates()
        return ax.plot_trisurf(*[xy[:,i] for i in range(gdim)],
                triangles=mesh.cells(), **kwargs)
    elif tdim == 1:
        x = [mesh.coordinates()[:,i] for i in range(gdim)]
        if gdim == 1:
            x.append(np.zeros_like(x[0]))
            ax.set_yticks([])
        marker = kwargs.pop('marker', 'o')
        return ax.plot(*x, marker=marker, **kwargs)
    else:
        assert False, "this code should not be reached"

# TODO: This is duplicated somewhere else
def create_cg1_function_space(mesh, sh):
    r = len(sh)
    if r == 0:
        V = dolfin.FunctionSpace(mesh, "CG", 1)
    elif r == 1:
        V = dolfin.VectorFunctionSpace(mesh, "CG", 1, dim=sh[0])
    else:
        V = dolfin.TensorFunctionSpace(mesh, "CG", 1, shape=sh)
    return V

def mplot_expression(ax, f, mesh, **kwargs):
    # TODO: Can probably avoid creating the function space here by restructuring
    #       mplot_function a bit so it can handle Expression natively
    V = create_cg1_function_space(mesh, f.ufl_shape)
    g = dolfin.interpolate(f, V)
    return mplot_function(ax, g, **kwargs)

def mplot_function(ax, f, **kwargs):
    mesh = f.function_space().mesh()
    gdim = mesh.geometry().dim()
    tdim = mesh.topology().dim()

    # Extract the function vector in a way that also works for subfunctions
    try:
        fvec = f.vector()
    except RuntimeError:
        fspace = f.function_space()
        try:
            fspace = fspace.collapse()
        # Happens for part of MultiMeshFunction; no way detecting elsewhere
        except RuntimeError:
            cpp.warning("Probably trying to plot MultiMeshFunction "
                        "part. Continuing without plotting...")
            return
        fvec = dolfin.interpolate(f, fspace).vector()

    if fvec.size() == mesh.num_cells():
        # DG0 cellwise function
        C = fvec.array() # NB! Assuming here dof ordering matching cell numbering
        if gdim == 2 and tdim == 2:
            return ax.tripcolor(mesh2triang(mesh), C, **kwargs)
        elif gdim == 3 and tdim == 2: # surface in 3d
            # FIXME: Not tested, probably broken
            xy = mesh.coordinates()
            shade = kwargs.pop("shade", True)
            return ax.plot_trisurf(mesh2triang(mesh), xy[:,2], C, shade=shade, **kwargs)
        elif gdim == 1 and tdim == 1:
            x = mesh.coordinates()[:,0]
            nv = len(x)
            # Insert duplicate points to get piecewise constant plot
            xp = np.zeros(2*nv-2)
            xp[0] = x[0]
            xp[-1] = x[-1]
            xp[1:2*nv-3:2] = x[1:-1]
            xp[2:2*nv-2:2] = x[1:-1]
            Cp = np.zeros(len(xp))
            Cp[0:len(Cp)-1:2] = C
            Cp[1:len(Cp):2] = C
            return ax.plot(xp, Cp, *kwargs)
        #elif tdim == 1: # FIXME: Plot embedded line
        else:
            raise AttributeError('Matplotlib plotting backend only supports 2D mesh for scalar functions.')

    elif f.value_rank() == 0:
        # Scalar function, interpolated to vertices
        # TODO: Handle DG1?
        C = f.compute_vertex_values(mesh)
        if gdim == 2 and tdim == 2:
            mode = kwargs.pop("mode", "color")
            if mode == "color":
                shading = kwargs.pop("shading", "flat")
                return ax.tripcolor(mesh2triang(mesh), C, shading=shading, **kwargs)
            elif mode == "warp":
                from matplotlib import cm
                cmap = kwargs.pop("cmap", cm.jet)
                linewidths = kwargs.pop("linewidths", 0)
                return ax.plot_trisurf(mesh2triang(mesh), C, cmap=cmap,
                        linewidths=linewidths, **kwargs)
            elif mode == "wireframe":
                return ax.triplot(mesh2triang(mesh), **kwargs)
            elif mode == "contour":
                return ax.tricontour(mesh2triang(mesh), C, **kwargs)
        elif gdim == 3 and tdim == 2: # surface in 3d
            # FIXME: Not tested
            from matplotlib import cm
            cmap = kwargs.pop("cmap", cm.jet)
            return ax.plot_trisurf(mesh2triang(mesh), C, cmap=cmap, **kwargs)
        elif gdim == 3 and tdim == 3:
            # Volume
            # TODO: Isosurfaces?
            # Vertex point cloud
            X = [mesh.coordinates()[:, i] for i in range(gdim)]
            return ax.scatter(*X, c=C, **kwargs)
        elif gdim == 1 and tdim == 1:
            x = mesh.coordinates()[:,0]
            ax.set_aspect('auto')

            # Setting limits for Line2D objects
            vmin = kwargs.pop("vmin", None)
            vmax = kwargs.pop("vmax", None)
            ax.set_ylim([vmin, vmax])
            return ax.plot(x, C, **kwargs)
        #elif tdim == 1: # FIXME: Plot embedded line
        else:
            raise AttributeError('Matplotlib plotting backend only supports 2D mesh for scalar functions.')

    elif f.value_rank() == 1:
        # Vector function, interpolated to vertices
        w0 = f.compute_vertex_values(mesh)
        nv = mesh.num_vertices()
        if len(w0) != gdim*nv:
            raise AttributeError('Vector length must match geometric dimension.')
        X = mesh.coordinates()
        X = [X[:, i] for i in range(gdim)]
        U = [w0[i*nv: (i+1)*nv] for i in range(gdim)]

        # Compute magnitude
        C = U[0]**2
        for i in range(1,gdim):
            C += U[i]**2
        C = np.sqrt(C)

        mode = kwargs.pop("mode", "glyphs")
        if mode == "glyphs":
            args = X + U + [C]
            if gdim == 3:
                # 3d quiver plot works only since matplotlib 1.4
                import matplotlib
                if StrictVersion(matplotlib.__version__) < '1.4':
                    cpp.warning('Matplotlib version %s does not support 3d '
                                'quiver plot. Continuing without plotting...'
                                % matplotlib.__version__)
                    return

                length = kwargs.pop("length", 0.1)
                return ax.quiver(*args, length=length, **kwargs)
            else:
                return ax.quiver(*args, **kwargs)
        elif mode == "displacement":
            Xdef = [X[i] + U[i] for i in range(gdim)]
            import matplotlib.tri as tri
            if gdim == 2 and tdim == 2:
                # FIXME: Not tested
                triang = tri.Triangulation(Xdef[0], Xdef[1], mesh.cells())
                shading = kwargs.pop("shading", "flat")
                return ax.tripcolor(triang, C, shading=shading, **kwargs)
            else:
                # Return gracefully to make regression test pass without vtk
                cpp.warning('Matplotlib plotting backend does not support '
                            'displacement for %d in %d. Continuing without '
                            'plotting...' % (tdim, gdim))
                return

def mplot_meshfunction(ax, obj, **kwargs):
    mesh = obj.mesh()
    tdim = mesh.topology().dim()
    d = obj.dim()
    if tdim == 2 and d == 2:
        C = obj.array()
        triang = mesh2triang(mesh)
        assert not kwargs.pop("facecolors", None), "Not expecting 'facecolors' in kwargs"
        return ax.tripcolor(triang, facecolors=C, **kwargs)
    #elif tdim == 1 and d == 1:
    #elif d == 0: # vertices
    else:
        # Return gracefully to make regression test pass without vtk
        cpp.warning('Matplotlib plotting backend does not support mesh '
                    'function of dim %d. Continuing without plotting...' % d)
        return

def mplot_dirichletbc(ax, obj, **kwargs):
    raise AttributeError("Matplotlib plotting backend doesn't handle DirichletBC.")

def _plot_matplotlib(obj, mesh, kwargs):
    if isinstance(obj, (cpp.MultiMesh, cpp.MultiMeshFunction, cpp.MultiMeshDirichletBC)):
        cpp.warning("Don't know how to plot type %s" % type(obj))
        return

    # Avoid importing pyplot until used
    try:
        import matplotlib.pyplot as plt
    except:
        print("*** Warning: matplotlib.pyplot not available, cannot plot")
        return

    gdim = mesh.geometry().dim()
    if gdim == 3 or kwargs.get("mode") in ("warp",):
        # Importing this toolkit has side effects enabling 3d support
        from mpl_toolkits.mplot3d import axes3d
        # Enabling the 3d toolbox requires some additional arguments
        ax = plt.gca(projection='3d')
    else:
        ax = plt.gca()
    ax.set_aspect('equal')

    title = kwargs.pop("title", None)
    if title is not None:
        ax.set_title(title)

    # Translate range_min/max kwargs supported by VTKPlotter
    vmin = kwargs.pop("range_min", None)
    vmax = kwargs.pop("range_max", None)
    if vmin and not "vmin" in kwargs:
        kwargs["vmin"] = vmin
    if vmax and not "vmax" in kwargs:
        kwargs["vmax"] = vmax

    # Let's stays consistent and use mode='color' instead of 'surface'
    if kwargs.get("mode") == "surface":
        cpp.deprecation("plot kwarg mode='surface'", "2016.1",
                        "Use mode='color' instead.")
        kwargs["mode"] = "color"

    # Drop unsupported kwargs and inform user
    _unsupported_kwargs = ["interactive", "rescale", "wireframe"]
    for kw in _unsupported_kwargs:
        if kwargs.pop(kw, None):
            cpp.warning("Matplotlib backend does not support '%s' kwarg yet. "
                        "Ignoring it..." % kw)

    if isinstance(obj, cpp.Function):
        return mplot_function(ax, obj, **kwargs)
    elif isinstance(obj, cpp.Expression):
        return mplot_expression(ax, obj, mesh, **kwargs)
    elif isinstance(obj, cpp.Mesh):
        return mplot_mesh(ax, obj, **kwargs)
    elif isinstance(obj, cpp.DirichletBC):
        return mplot_dirichletbc(ax, obj, **kwargs)
    elif isinstance(obj, _meshfunction_types):
        return mplot_meshfunction(ax, obj, **kwargs)
    else:
        raise AttributeError('Failed to plot %s' % type(obj))


# Compatibility with book
def _VTKPlotter_write_ps(self, *args, **kwargs) :
    print("*** Warning: VTKPlotter::write_ps() is not implemented -- use write_pdf instead")

_objects_referenced_from_plot_windows = {}
def _plot_cpp(object, mesh, kwargs):
    # Convert kwargs to cpp format
    p = cpp.Parameters()
    for key in kwargs:
        try:
            p.add(key, kwargs[key])
        except TypeError:
            cpp.warning("Incompatible type for keyword argument \"%s\". Ignoring." % key)

    if isinstance(object, cpp.Expression):
        plot_object = cpp.plot(object, mesh, p)
    elif isinstance(object, cpp.MultiMesh):
        return cpp.io._plot_multimesh(object)
    else:
        plot_object = cpp.plot(object, p)

    # Compatibility with book
    plot_object.write_ps = _VTKPlotter_write_ps

    # Avoid premature deletion of plotted objects if they go out of scope
    # before the plot window is closed. The plotter itself is safe, since it's
    # created in the plot() C++ function, not directly from Python. But the
    # Python plotter proxy may disappear, so we can't store the references
    # there.
    global _objects_referenced_from_plot_windows
    _objects_referenced_from_plot_windows[plot_object.key()] = (object, mesh, p)

    return plot_object

def plot(object, *args, **kwargs):
    """
    Plot given object.

    *Arguments*
        object
            a :py:class:`Mesh <dolfin.cpp.Mesh>`, a :py:class:`MeshFunction
            <dolfin.cpp.MeshFunction>`, a :py:class:`Function
            <dolfin.functions.function.Function>`, a :py:class:`Expression`
            <dolfin.cpp.Expression>, a :py:class:`DirichletBC`
            <dolfin.cpp.DirichletBC>, a :py:class:`FiniteElement
            <ufl.FiniteElement>`, or a :py:class:`MultiMesh
            <dolfin.cpp.MultiMesh>`.

    *Examples of usage*
        In the simplest case, to plot only e.g. a mesh, simply use

        .. code-block:: python

            mesh = UnitSquare(4,4)
            plot(mesh)

        Use the ``title`` argument to specify title of the plot

        .. code-block:: python

            plot(mesh, tite="Finite element mesh")

        It is also possible to plot an element

        .. code-block:: python

            element = FiniteElement("BDM", tetrahedron, 3)
            plot(element)

        Vector valued functions can be visualized with an alternative mode

        .. code-block:: python

            plot(u, mode = "glyphs")

        A more advanced example

        .. code-block:: python

            plot(u,
                 wireframe = True,              # use wireframe rendering
                 interactive = False,           # do not hold plot on screen
                 scalarbar = False,             # hide the color mapping bar
                 hardcopy_prefix = "myplot",    # default plotfile name
                 scale = 2.0,                   # scale the warping/glyphs
                 title = "Fancy plot",          # set your own title
                 backend = "vtk"                # choose plotting backend
                 )

    """

    # Plot element
    if isinstance(object, ufl.FiniteElementBase):
        if os.environ.get("DOLFIN_NOPLOT", "0") != "0":
            return
        import ffc
        return ffc.plot(object, *args, **kwargs)

    # Get mesh from explicit mesh kwarg, only positional arg, or via object
    mesh = kwargs.pop('mesh', None)
    if isinstance(object, cpp.Mesh):
        if mesh is not None and mesh.id() != object.id():
            cpp.dolfin_error("plotting.py",
                             "plot mesh",
                             "Got different mesh in plot object and keyword argument")
        mesh = object
    if mesh is None:
        if isinstance(object, cpp.Function):
            mesh = object.function_space().mesh()
        elif hasattr(object, "mesh"):
            mesh = object.mesh()

    # Expressions do not carry their own mesh
    if isinstance(object, cpp.Expression) and mesh is None:
        cpp.dolfin_error("plotting.py",
                         "plot expression",
                         "Expecting a mesh as keyword argument")

    # Try to project if object is not a standard plottable type
    if not isinstance(object, _plottable_types):
        from dolfin.fem.projection import project
        try:
            cpp.info("Object cannot be plotted directly, projecting to"\
                     " piecewise linears.")
            object = project(object, mesh=mesh)
            mesh = object.function_space().mesh()
        except Exception as e:
            msg = "Don't know how to plot given object:\n  %s\n" \
                  "and projection failed:\n  %s" % (str(object), str(e))
            cpp.dolfin_error("plotting.py", "plot object", msg)

    return _plot_matplotlib(object, mesh, kwargs)
    # Select backend
    '''
    backend = kwargs.pop("backend", None) or cpp.parameters["plotting_backend"]
    if backend == "vtk":
        return _plot_cpp(object, mesh, kwargs)
    elif backend == "matplotlib":
        return _plot_matplotlib(object, mesh, kwargs)
    elif backend == "none":
        cpp.warning("Plotting backend set to 'none'. Plotting disabled.")
        return
    else:
        cpp.dolfin_error("plotting.py", "plot object",
                         "Unknown plotting backend '%s'" % backend)
    '''