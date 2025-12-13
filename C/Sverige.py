from dolfin import *
from numpy import dtype

# Create mesh and define function space
mesh = Mesh("C/meshes/sweden_coarse.xml.gz")
mesh = Mesh("C/meshes/sweden_coarse.xml.gz")
import numpy as np

coords = mesh.coordinates()

xmin = coords[:, 0].min()
xmax = coords[:, 0].max()

ymin = coords[:, 1].min()
ymax = coords[:, 1].max()

print(xmin, xmax)
print(ymin, ymax)
# # Construct the finite element space
# P1 = FiniteElement("Lagrange", mesh.ufl_cell() , 1)
# TH = P1 * P1 * P1
# W = FunctionSpace(mesh , TH)

# # Define parameters :
# T = 500
# dt = 0.5
# t = 0
# delta1 = 1
# delta2 = 1
# delta3 = 1
# alpha = 0.4
# beta = 0.8
# gamma = 0.8
# zeta = 2
# L_0 = 0.4
# l = 0.6
# m = 0.12

# # Class representing the intial conditions
# class InitialConditions(UserExpression):
#     def eval(self, values , x) :
#         values[1] = 4/15 - 2*10**(-7)*(x[0]-0.1*x[1]-350)**2
#         values[0] = 0
#         values[2] = 22/45 - 3*10**(-5) * (x[0]-450) - 1.2*10**(-4) * (x[1]-15)
#     def value_shape(self) :
#         return (3 ,)

# # Define the nonlinear term
# def N(u_h):
#     u = u_h[0]
#     v = u_h[1]
#     w = u_h[2]
#     N1 = - alpha * u * (1 - u/(L_0 + l*v))
#     N2 = - beta * v * (1 - v) + v*w / (alpha + v + m * u)
#     N3 = gamma * w - zeta * v * w / (alpha + v + m * u)
#     return as_vector([N1 , N2 , N3])

# # Define initial condition
# indata = InitialConditions(degree = 2)
# u0 = Function(W)
# u0 = interpolate(indata , W)

# # Create bilinear and linear forms
# psi = TestFunction(W)
# u_old = Function(W)
# u_new = Function(W)
# u = TrialFunction(W)
# D = as_matrix([[delta1 , 0 , 0] ,
#                    [0 , delta2 , 0] ,
#                    [0 , 0 , delta3]])

# # Bilinear form: crank nicolson
# F = inner((u - u_old)/dt, psi)*dx + inner(N(u_old), psi)*dx +inner(D*grad((u+u_old)/2), grad(psi))*dx
# a = lhs(F)
# L = rhs(F)

# # Set an output file
# file = File("C/Solutions/sverige_solution.pvd")

# # Set initial condition
# u_old.assign(u0)

# # Open CSV file for population data
# csv_file = open("C/Solutions/sverige_populations.csv", "w")
# csv_file.write("time,population_u,population_v,population_w\n")

# # Time - stepping
# while t < T:
#     # Define the integrals
#     M0 = u_old[0] * dx
#     M1 = u_old[1] * dx
#     M2 = u_old[2] * dx
#     # compute the functional
#     population_u = assemble ( M0 )
#     population_v = assemble ( M1 )
#     population_w = assemble ( M2 )
    
#     # Write to CSV
#     csv_file.write(f"{t},{population_u},{population_v},{population_w}\n")
    
#     if int(t) % 100 == 0:
#         print("Time : ", t)
#         file << u_old

#     solve(a == L, u_new)
#     u_old.assign(u_new)
#     t += dt

# csv_file.close()
