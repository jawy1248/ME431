#%%

# Imports
import sympy as sp
from sympy import eye, sin, cos, diff, Matrix, symbols, Function, pretty_print, simplify, init_printing, latex
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
from IPython.display import Math, display

# Define symbols
t, m1, m2, ell, g = symbols('t, m1, m2, ell, g')
z = dynamicsymbols('z')
theta = dynamicsymbols('theta')

# Define generalized coords
q = Matrix([[z], [theta]])
qdot = q.diff(t)

# Get position and velocity
p1 = Matrix([[z*cos(theta)], [z*sin(theta)], [0]])
p2 = Matrix([[0.5*ell*cos(theta)], [0.5*ell*sin(theta)], [0]])
v1 = p1.diff(t)
v2 = p2.diff(t)

# Define J and omega
omega = Matrix([[0], [0], [theta.diff(t)]])
J = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, m2*ell**2/12.0]])

# Get the kinetic energy
K = simplify(0.5*m1*v1.T*v1 + 0.5*m2*v2.T*v2 + 0.5*omega.T*J*omega)
K = K[0,0]
display("KE:")
display(Math(vlatex(K)))

# Get the potential energy
P = m1*g*z*sin(theta) + m2*g*0.5*ell*sin(theta)

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
thetad = theta.diff(t)
thetadd = thetad.diff(t)
F = symbols('F')

# Get the RHS
RHS = Matrix([[0], [F*ell*cos(theta)]])

# Get whole EOM
FullEOM = LHS - RHS

# Solve and display result
result = simplify(sp.solve(FullEOM, (zdd, thetadd)))
zdd_eom = result[zdd]
thetadd_eom = result[thetadd]

eoms = Matrix([[zdd_eom],[thetadd_eom]])
display("EOMs: ")
display(Math(vlatex(simplify(eoms))))

gTerm = (3*g*m1*z*cos(theta))/(ell**2*m2 + 3*m1*z**2)
svf = Matrix([[zdd_eom], [simplify(thetadd_eom + gTerm)], [zd], [thetad]])
states = Matrix([[zd], [thetad], [z], [theta]])
inputs = Matrix([[F]])

A = svf.jacobian(states)
B = svf.jacobian(inputs)

A_lin = simplify(A.subs([(thetad, 0.0), (zd, 0.0), (theta, 0), (F, (2*m1*g*z + ell*m2*g)/(2*ell))]))
B_lin = simplify(B.subs([(thetad, 0.0), (zd, 0.0), (theta, 0), (F, (2*m1*g*z + ell*m2*g)/(2*ell))]))

A_full = simplify(A.subs([(thetad, 0.0), (zd, 0.0), (theta, 0), (F, (2*m1*g*z + ell*m2*g)/(2*ell)), (m1, 0.35), (m2, 2), (ell, 0.5), (g, 9.8), (z, 0)]))
B_full = simplify(B.subs([(thetad, 0.0), (zd, 0.0), (theta, 0), (F, (2*m1*g*z + ell*m2*g)/(2*ell)), (m1, 0.35), (m2, 2), (ell, 0.5), (g, 9.8), (z, 0)]))

display("Linear EOMs (A) then (B):")
display(Math(vlatex(A_full)))
display(Math(vlatex(B_full)))

# Getting T.F.
C = Matrix([[0, 0, 1, 0], [0, 0, 0, 1]])
D = Matrix([[0], [0]])
I = eye(4)
s = symbols('s')

TF = simplify(C @ (s*I - A_lin).inv() @ B_lin + D)

display("TF:")
TF = simplify(TF.subs([(z, 0.0)]))
display(Math(vlatex(simplify(TF))))
display(Math(vlatex(simplify(TF[0]/TF[1]))))

# %%
