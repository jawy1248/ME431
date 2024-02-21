#%%

# Imports
import sympy as sp
from sympy import eye, sin, cos, diff, Matrix, symbols, Function, pretty_print, simplify, init_printing, latex
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
eoms = Matrix([[zdd_eom]])
display("EOMs:")
display(Math(vlatex(eoms)))

svf = Matrix([[zd],[zdd_eom]])
states = Matrix([[z], [zd]])
inputs = Matrix([[F]])

A = svf.jacobian(states)
B = svf.jacobian(inputs)

A_lin = simplify(A.subs([(zd, 0.0), (F, z*k), (z, 0), (m, 5), (k, 3), (b, 0.5)]))
B_lin = simplify(B.subs([(zd, 0.0), (F, z*k), (z, 0), (m, 5), (k, 3), (b, 0.5)]))

display("Linear EOMs (A) then (B):")
display(Math(vlatex(A_lin)))
display(Math(vlatex(B_lin)))

# Getting T.F.
C = Matrix([[1, 0]])
D = Matrix([[0]])
I = eye(2)
s = symbols('s')

TF = simplify(C @ (s*I - A).inv() @ B + D)

display("TF:")
display(Math(vlatex(TF)))

# %% 
# Solve the block diagram
k_p, k_d = symbols('k_p, k_d')
# charEQ = 1/TF[0]
charEQ = s**2 + (0.1 + k_d/5)*s + (0.6 + k_p/5)
sVar = sp.solve(charEQ, s)

p1 = -1.0
p2 = -1.5

eq1 = sVar[1] - p1
eq2 = sVar[0] - p2

eqn = Matrix([[eq1], [eq2]])

res = sp.solve(eqn, k_d, k_p)
display(Math(vlatex(res)))

# %%
