from dolfin import *

# Create mesh and define function space
mesh = Mesh("C/meshes/circle.xml.gz")

# Construct the finite element space
P1 = FiniteElement("Lagrange", mesh.ufl_cell() , 1)
TH = P1 * P1 * P1
W = FunctionSpace(mesh , TH)

# Define parameters :
T = 500
dt = 0.5
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
        values[0] = 0
        values[1] = 4/15 - 2*10**(-7)*(x[0]-0.1*x[1]-350)**2
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
u = Function(W)
D = as_matrix([[delta1 , 0 , 0] ,
                   [0 , delta2 , 0] ,
                   [0 , 0 , delta3]])
a = dot(u, psi)*dx + dot(N(u), psi)*dx + inner(D*grad(u), grad(psi))*dx
L = 0

# Set an output file
file = File("C/Solutions/solution.pvd")

# Compute solution
solve(a == L, u,
    solver_parameters={"newton_solver": {"relative_tolerance": 1e-6}})

file << u
# # Set initial condition
# ...

# # Time - stepping
# while t < T :
#     # assign u0
#     u0.assign(u)
#     ...
#     ...