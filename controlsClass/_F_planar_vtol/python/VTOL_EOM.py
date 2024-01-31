#%%

# Imports
import sympy as sp
from sympy import sin, cos, diff, Matrix, symbols, Function, pretty_print, simplify, init_printing, latex
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
from IPython.display import Math, display

# Define symbols
t, m1, m2, d, g, J1, mu = symbols('t, m1, m2, d, g, J1, mu')
z = dynamicsymbols('z')
h = dynamicsymbols('h')
theta = dynamicsymbols('theta')

# Define generalized coords
q = Matrix([[z], [h], [theta]])
qdot = q.diff(t)

# Get position and velocity
pl = Matrix([[z - d*cos(theta)], [h - d*sin(theta)], [0]])
pc = Matrix([[z], [h], [0]])
pr = Matrix([[z + d*cos(theta)], [h + d*sin(theta)], [0]])

vl = pl.diff(t)
vc = pc.diff(t)
vr = pr.diff(t)

# Define J and omega
omega = Matrix([[0], [0], [theta.diff(t)]])
J = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, J1]])

# Get the kinetic energy
K = simplify(0.5*m1*vl.T*vl + 0.5*m2*vc.T*vc + 0.5*m1*vr.T*vr + 0.5*omega.T*J*omega)
K = K[0,0]
display("KE:")
display(Math(vlatex(K)))

#%%
# Get the potential energy
P = (m1 + m2 + m1)*g*h

# Get the Lagrangian
L = simplify(K-P)

#  Get the equation for the LHS
LHS = simplify((L.diff(qdot)).diff(t) - L.diff(q))
display("LHS:")
display(Math(vlatex(LHS)))

#%%
# Getting generalized forces and damping
# Define more variables
zd = z.diff(t)
zdd = zd.diff(t)
hd = h.diff(t)
hdd = hd.diff(t)
thetad = theta.diff(t)
thetadd = thetad.diff(t)
f1, f3 = symbols('f1, f3')

# Get the RHS
RHS = Matrix([[-sin(theta)*(f1+f3) - mu*zd], [cos(theta)*(f1+f3)], [d*(f3-f1)]])

# Get whole EOM
FullEOM = LHS - RHS

# Solve and display result
result = simplify(sp.solve(FullEOM, (zdd, hdd, thetadd)))
zdd_eom = result[zdd]
hdd_eom = result[hdd]
thetadd_eom = result[thetadd]

eoms = Matrix([[zdd_eom],[hdd_eom],[thetadd_eom]])
display(Math(vlatex(eoms)))

# %%
