__author__ = 'svetlana'

from dolfin import *
from dolfin.cpp.common import tic, toc, set_log_level
from dolfin.cpp.mesh import Mesh, SubDomain, FacetFunction

import math
import os
import scipy.io
import scipy.optimize
import scipy.stats
import random
import numpy

# Adjust log level
set_log_level(PROGRESS)

# Turn on optimization
parameters["form_compiler"]["cpp_optimize"] = True

def get_list_of_power_type_basis_functions(N):
    phi_list = []
    for jj in xrange(0, N+1):
        if jj == 0:
            for ii in xrange(1, N+1):
                phi_list.append(Expression("- 1.0 / (i + 1) + pow(x[0], i)", i=ii))
        else:
            for ii in xrange(0, N+1):
                phi_list.append(Expression("pow(x[0], i) * pow(x[1], j)", i=ii, j=jj))
    return phi_list

def get_list_of_fourier_type_basis_functions(N):
    phi_list = []
    for ii in xrange(0, N+1):
        if ii == 0:
            for jj in xrange(1, N+1):
                phi_list.append(Expression("- 1.0 + cos(pi*j*x[1])", j=jj))
        else:
            for jj in xrange(0, N+1):
                phi_list.append(Expression("cos(pi*i*x[0]) * cos(pi*j*x[1])", i=ii, j=jj))
    return phi_list

# Import functionality from the folder 'lib'
import os, sys

lib_path = os.path.abspath('lib')
sys.path.append(lib_path)
import postprocess

# ---------------------------------------------------------------------------------#
# Calculate reference triangle constant C_hat_G and C_hat_P

tic()
h = 1

def f1_tan(x):
    return x*cos(x)/sin(x) + 1


def f2_tan(x):
    return tan(x) + tanh(x)

zeta_1 = scipy.optimize.newton(f1_tan, DOLFIN_PI/2)
zeta_2 = scipy.optimize.newton(f2_tan, DOLFIN_PI)

sigma_1 = abs(zeta_2 * tan(zeta_2))

C_hat_P = h / zeta_1
C_hat_G = sqrt(h / sigma_1)

# Calculate C_hat according to the formula
#rho = [0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]
rho = [1/1]
alpha = DOLFIN_PI / 2
mu = [0.5 * (1 + x ** 2 + math.sqrt(1 + x ** 4 + 2 * cos(2 * alpha) * x ** 2)) for x in rho]
c_rho_alpha = []
c_rho = []

C_G = []
C_P = []

for i in range(0, len(rho)):
    c_rho_alpha.append(sqrt(mu[i] / (sin(alpha) * rho[i])))
    c_rho.append(sqrt(mu[i]))
    C_G.append(c_rho_alpha[i] * C_hat_G * sqrt(h))
    C_P.append(c_rho[i] * C_hat_P * h)

print "rho = ", rho
print "c(rho, alpha) = ", c_rho_alpha
print "C_hat_G = ", C_hat_G
print "C_G = ", C_G

print "c(rho) = ", c_rho
print "C_hat_P = ", C_hat_P
print "C_P = ", C_P

if alpha == DOLFIN_PI / 2:
    # Define the name of xml-file
    #file_names = ['alpha-pi2-h-14-v-145.xml', 'alpha-pi2-h-12-v-289.xml', 'alpha-pi2-h-1-v-289.xml', 'alpha-pi2-h-32-v-181.xml', 'alpha-pi2-h-2-v-217.xml']
    file_names = ['alpha-pi2-h-1-v-289.xml']
    # Define the file path
    file_path = '/data/pi-2/'

elif alpha == DOLFIN_PI / 4:
    # Define the name of xml-file
    #file_names = ['alpha-pi4-h-14.xml', 'alpha-pi4-h-12.xml', 'alpha-pi4-h-1.xml', 'alpha-pi4-h-32.xml', 'alpha-pi4-h-2.xml']
    file_names = ['alpha-pi4-h-1.xml']
    # Define the file path
    file_path = '/data/pi-4/'

elif alpha == DOLFIN_PI / 3:
    # Define the name of xml-file
    #file_names = ['alpha-pi3-h-14.xml', 'alpha-pi3-h-12.xml', 'alpha-pi3-h-1.xml', 'alpha-pi3-h-32.xml', 'alpha-pi3-h-2.xml']
    file_names = ['alpha-pi3-h-1.xml']
    # Define the file path
    file_path = '/data/pi-3/'

elif alpha == DOLFIN_PI / 6:
    # Define the name of xml-file
    #file_names = ['alpha-pi6-h-sqrt2-2.xml', 'alpha-pi6-h-12.xml', 'alpha-pi6-h-1.xml', 'alpha-pi6-h-32.xml', 'alpha-pi6-h-2.xml']
    file_names = ['alpha-pi6-h-1.xml']
    # Define the file path
    file_path = '/data/pi-6/'

elif alpha == DOLFIN_PI / 10:

    file_names = ['alpha-pi10-h-1.xml']
    # Define the file path
    file_path = '/data/pi-10/'

elif alpha == DOLFIN_PI / 30:

    file_names = ['alpha-pi30-h-1.xml']
    # Define the file path
    file_path = '/data/pi-30/'

elif alpha == 3 * DOLFIN_PI / 4:

    file_names = ['alpha-3pi4-h-1.xml']
    # Define the file path
    file_path = '/data/3-pi-4/'

elif alpha == 2 * DOLFIN_PI / 3:

    file_names = ['alpha-2pi3-h-1.xml']
    # Define the file path
    file_path = '/data/2-pi-3/'

elif alpha == 5 * DOLFIN_PI / 6:

    file_names = ['alpha-5pi6-h-1.xml']
    # Define the file path
    file_path = '/data/5-pi-6/'

elif alpha == 9 * DOLFIN_PI / 10:

    file_names = ['alpha-9pi10-h-1.xml']
    # Define the file path
    file_path = '/data/9-pi-10/'

elif alpha == 29 * DOLFIN_PI / 30:

    file_names = ['alpha-29pi30-h-1.xml']
    # Define the file path
    file_path = '/data/29-pi-30/'

# Maximal N in the power series or in the fourier series
N_power = numpy.array([6, 7])

for file_name in file_names: # if the there are list of xml-file we run through the loop

    # Load mesh from the xml-file
    project_path = postprocess.get_project_path()
    mesh = Mesh(project_path +
                file_path +
                file_name)
    print "mesh: %s, h_max = %f" % (file_name, mesh.hmax())


    class ZeroBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and x[1] <= DOLFIN_EPS

    # Mark boundaries
    zero_boundary = ZeroBoundary()
    boundaries = FacetFunction("uint", mesh)
    boundaries.set_all(0)
    zero_boundary.mark(boundaries, 1)

    # Redefine boundary measure
    ds = Measure('ds')[boundaries]

    V = FunctionSpace(mesh, "Lagrange", 3)

    for n_power in range(N_power.size):

        N = N_power[n_power]
        phi = get_list_of_power_type_basis_functions(N)
        #phi = get_list_of_fourier_type_basis_functions(N)

        M = len(phi)

        K = numpy.array([[0.0 for x in xrange(M)] for x in xrange(M)])
        S = numpy.array([[0.0 for x in xrange(M)] for x in xrange(M)])
        T = numpy.array([[0.0 for x in xrange(M)] for x in xrange(M)])

        for j in range(0, M):
            for i in range(j, M):
                phi_i_ = phi[i]
                phi_j_ = phi[j]

                phi_i = interpolate(phi_i_, V)
                phi_j = interpolate(phi_j_, V)

                k = assemble(phi_j * phi_i * ds(1))
                t = assemble(phi_j * phi_i * dx)
                s = assemble(inner(grad(phi_j), grad(phi_i)) * dx)

                K[j, i] = k
                T[j, i] = t
                S[j, i] = s

                # Symmetric matrix
                if i > j:
                    K[i, j] = k
                    S[i, j] = s
                    T[i, j] = t
                #print "iter %d, %d: %f, %f" % (j, i, K[j, i], S[j, i])

        eigs_cp, v_cp = scipy.linalg.eig(T, S)
        eigs_cg, v_cg = scipy.linalg.eig(K, S)

        idx_cp = eigs_cp.argsort()
        sorted_eigs_cp = eigs_cp[idx_cp]
        sorted_v_cp = v_cp[:, idx_cp]

        idx_cg = eigs_cg.argsort()
        sorted_eigs_cg = eigs_cp[idx_cg]
        sorted_v_cg = v_cg[:, idx_cg]


        cp_1 = 1/sqrt(abs(sorted_eigs_cp[0]))
        cp_2 = 1/sqrt(abs(sorted_eigs_cp[1]))
        v_cp_1 = sorted_v_cp[:, 0]
        v_cp_2 = sorted_v_cp[:, 1]

        cg_1 = 1/sqrt(abs(sorted_eigs_cg[0]))
        cg_2 = 1/sqrt(abs(sorted_eigs_cg[1]))
        v_cg_1 = sorted_v_cg[:, 0]
        v_cg_2 = sorted_v_cg[:, 1]



    print 'N %d, M %d:' % (N, M)
    print 'cp: %e > %e' % (cp_1, cp_2)
    print 'cg: %e > %e' % (cg_1, cg_2)

path = postprocess.construct_result_path('FEniCS-functions-power-alpha-' + str(alpha) + '-rho-' + str(rho[0]) + '-N-' + str(N))
scipy.io.savemat(project_path + path, {"cp_1": cp_1, "cp_2": cp_2, "cg_1": cg_1, "cg_2": cg_2, "v_cp_1": v_cp_1, "v_cp_2": v_cp_2, "v_cg_1": v_cg_1, "v_cg_2": v_cg_2})
