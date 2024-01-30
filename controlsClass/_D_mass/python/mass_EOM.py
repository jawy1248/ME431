#%%

# Imports
import sympy as sp
from sympy import sin, cos, diff, Matrix, symbols, Function, pretty_print, simplify, init_printing, latex
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
from IPython.display import Math, display

# Define symbols
t, m, k, b = symbols('t, m, k, b')
z = dynamicsymbols('z')

# Define generalized coords
q = Matrix([[z]])
qdot = q.diff(t)

# Get position and velocity
p = Matrix([[z], [0], [0]])
v = p.diff(t)

# Get the kinetic energy
K = simplify(0.5*m*v.T*v)
K = K[0,0]

# Get the potential energy
P = 0.5*k*(z**2)

# Get the Lagrangian
L = simplify(K-P)

#  Get the equation for the LHS
LHS = simplify((L.diff(qdot)).diff(t) - L.diff(q))
display("KE:")
display(Math(vlatex(LHS)))

# Getting generalized forces and damping
# Define more variables
zd = z.diff(t)
zdd = zd.diff(t)
F, b = symbols('F, b')

# Get the RHS
RHS = Matrix([[F - b*zd]])

# Get whole EOM
FullEOM = LHS - RHS

# Solve and display result
result = simplify(sp.solve(FullEOM, (zdd)))
zdd_eom = result[zdd]
display("EOM:")
display(Math(vlatex(zdd_eom)))
# %%
