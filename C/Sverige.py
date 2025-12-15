from dolfin import *
import random
set_log_level(LogLevel.WARNING)

# Create mesh and define function space
mesh = Mesh("C/meshes/sweden.xml.gz")

# Construct the finite element space
P1 = FiniteElement("Lagrange", mesh.ufl_cell() , 1)
TH = P1 * P1 * P1
W = FunctionSpace(mesh , TH)

# Define parameters :
T = 1200
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

coords = mesh.coordinates()
xmin = coords[:, 0].min()
xmax = coords[:, 0].max()
ymin = coords[:, 1].min()
ymax = coords[:, 1].max()
print('Xmin and xmax: ',xmin, xmax)
print('Ymin and ymax: ', ymin, ymax)
Sundsvall_y = ymax - (ymax - ymin)/2
print("Sundsvall y-coordinate: ", Sundsvall_y)

def initial_values(species,x , my_number):
    if x < my_number:
        if species == 'u':
            return (5/1000) * random.uniform(0, 1)
        elif species == 'v':
            return 0.5 *(1 - random.uniform(0, 1))
        elif species == 'w':
            return 1/4 + (1/2) * random.uniform(0, 1)
        else:
            raise ValueError("Unknown species")
    else:
        return 1/100

# Class representing the intial conditions
class InitialConditions(UserExpression):
    def eval(self, values , x) :
        values[0] = initial_values('u', x[1], Sundsvall_y)
        values[1] = initial_values('v', x[1], Sundsvall_y)
        values[2] = initial_values('w', x[1], Sundsvall_y)
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

# Define the integrals
M0 = u_old[0] * dx
M1 = u_old[1] * dx
M2 = u_old[2] * dx

# Set an output file
name = 'sverige'
file = File(f"C/Solutions/{name}/{name}_solution.pvd")
# Open CSV file for population data
csv_file = open(f"C/Solutions/{name}/{name}_population.csv", "w")
csv_file.write("time,population_u,population_v,population_w\n")

# Set initial condition
u_old.assign(u0)
# Time - stepping
while t < T:
    print("Current time: ", t, end='\r')
    # compute the functional
    population_u = assemble ( M0 )
    population_v = assemble ( M1 )
    population_w = assemble ( M2 )
    
    # Write to CSV
    csv_file.write(f"{t},{population_u},{population_v},{population_w}\n")
    
    if t % 50 == 0:
        print("Saving solution at time: ", t)
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
