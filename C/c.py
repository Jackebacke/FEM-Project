from fileinput import filename
import os
from dolfin import *

# Create mesh and define function space
mesh = Mesh("C/meshes/circle.xml.gz")

# Construct the finite element space
P1 = FiniteElement("Lagrange", mesh.ufl_cell() , 1)
TH = P1 * P1 * P1
W = FunctionSpace(mesh , TH)

# Define parameters :
T = 1000
dt = 0.5
t = 0
delta1 = 1
delta2 = 1
delta3 = 1
alpha = 0.4
beta = 0.8
gamma = 0.8
zeta = 2
L_0 = 0.4
l = 0.6
m = 0.12

# Class representing the intial conditions
class InitialConditions(UserExpression):
    def eval(self, values , x) :
        values[1] = 4/15 - 2*10**(-7)*(x[0]-0.1*x[1]-350)**2
        values[0] = 0 # Part C1a)
        # values[0] = 0.1 * values[1] # Part C1b)
        values[2] = 22/45 - 3*10**(-5) * (x[0]-450) - 1.2*10**(-4) * (x[1]-15)
    def value_shape(self) :
        return (3 ,)

# Define the nonlinear term
def N(u_h):
    u = u_h[0]
    v = u_h[1]
    w = u_h[2]
    N1 = - alpha * u * (1 - u/(L_0 + l*v))
    N2 = - beta * v * (1 - v) + v*w / (alpha + v + m * u)
    N3 = gamma * w - zeta * v * w / (alpha + v + m * u)
    return as_vector([N1 , N2 , N3])

# Define initial condition
indata = InitialConditions(degree = 2)
u0 = Function(W)
u0 = interpolate(indata , W)

# Create bilinear and linear forms
psi = TestFunction(W)
u_old = Function(W)
u_new = Function(W)
u = TrialFunction(W)
D = as_matrix([[delta1 , 0 , 0] ,
                   [0 , delta2 , 0] ,
                   [0 , 0 , delta3]])
# Bilinear form: crank nicolson
F = inner((u - u_old)/dt, psi)*dx + inner(N(u_old), psi)*dx +inner(D*grad((u+u_old)/2), grad(psi))*dx
a = lhs(F)
L = rhs(F)

# Define integral forms once
M0 = u_old[0] * dx
M1 = u_old[1] * dx
M2 = u_old[2] * dx

# Output files
path = "C/Solutions/c1a"
file_name = "c1a"
file = File(os.path.join(path, f"{file_name}.pvd"))
# Open CSV file for population data
csv_file = open(os.path.join(path, f"{file_name}_population.csv"), "w")
csv_file.write("time,population_u,population_v,population_w\n")


# Set initial condition
u_old.assign(u0)
# Time - stepping
while t < T:
    # compute the functional (u_old changes each iteration)
    population_u = assemble(M0)
    population_v = assemble(M1)
    population_w = assemble(M2)
    
    # Write to CSV
    csv_file.write(f"{t},{population_u},{population_v},{population_w}\n")
    
    if t % 100 == 0:  # Save when t crosses multiples of 100
        print("Time : ", t)
        file << u_old

    solve(a == L, u_new)
    u_old.assign(u_new)
    t += dt

# compute the functional at final time step
population_u = assemble(M0)
population_v = assemble(M1)
population_w = assemble(M2)
# Write to CSV final time step
csv_file.write(f"{t},{population_u},{population_v},{population_w}\n")

csv_file.close()

# Save final solution
print(f"Final time: {t}")
file << u_old
