import numpy as np
import matplotlib.pyplot as plt
from dolfin import *

def initial_conditions(x, values):
    values[0] = 0
    values[1] = 4/15 - 2*10**(-7)*(x[0]-0.1*x[1]-350)**2
    values[2] = 22/45 - 3*10**(-5) * (x[0]-450) - 1.2*10**(-4) * (x[1]-15)

# Load the mesh from c.py
mesh = Mesh("C/meshes/circle.xml.gz")

# Create function space and interpolate initial conditions
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = P1 * P1 * P1
W = FunctionSpace(mesh, TH)

# Define and interpolate initial conditions
class InitialConditions(UserExpression):
    def eval(self, values, x):
        values[0] = 0
        values[1] = 4/15 - 2*10**(-7)*(x[0]-0.1*x[1]-350)**2
        values[2] = 22/45 - 3*10**(-5) * (x[0]-450) - 1.2*10**(-4) * (x[1]-15)
    def value_shape(self):
        return (3,)

indata = InitialConditions(degree=2)
u0 = interpolate(indata, W)

# Extract components using sub-function spaces
V0 = W.sub(0)
V1 = W.sub(1)

# Create functions for each component
u0_0 = Function(V0)
u0_1 = Function(V1)

# Assign the components
assign(u0_0, u0.sub(0))
assign(u0_1, u0.sub(1))

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 4))

# Plot first component
im0 = plot(u0_0, ax=axes[0], cmap='viridis')
axes[0].set_title('Component 1 (Initial Condition)')

# Plot second component
im1 = plot(u0_1, ax=axes[1], cmap='viridis')
axes[1].set_title('Component 2 (Initial Condition)')

plt.tight_layout()
plt.savefig('C/initial_conditions.png', dpi=150)
plt.show()
        