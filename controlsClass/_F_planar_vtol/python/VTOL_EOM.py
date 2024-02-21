#%%

# Imports
import sympy as sp
from sympy import eye, sin, cos, diff, Matrix, symbols, Function, pretty_print, simplify, init_printing, latex, monic
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

# Get the potential energy
P = (m1 + m2 + m1)*g*h

# Get the Lagrangian
L = simplify(K-P)

#  Get the equation for the LHS
LHS = simplify((L.diff(qdot)).diff(t) - L.diff(q))
display("LHS:")
display(Math(vlatex(LHS)))

# Getting generalized forces and damping
# Define more variables
zd = z.diff(t)
zdd = zd.diff(t)
hd = h.diff(t)
hdd = hd.diff(t)
thetad = theta.diff(t)
thetadd = thetad.diff(t)
# f1, f3 = symbols('f1, f3')
F, tau = symbols('F, tau')

# Get the RHS
RHS = Matrix([[-sin(theta)*F - mu*zd], [cos(theta)*F], [tau]])

# Get whole EOM
FullEOM = LHS - RHS

# Solve and display result
result = simplify(sp.solve(FullEOM, (zdd, hdd, thetadd)))
zdd_eom = result[zdd]
hdd_eom = result[hdd]
thetadd_eom = result[thetadd]

eoms = Matrix([[zdd_eom],[hdd_eom],[thetadd_eom]])
display(Math(vlatex(eoms)))

# Full
svf = Matrix([[zdd_eom], [hdd_eom], [thetadd_eom], [zd], [hd], [thetad]])
states = Matrix([[zd], [hd], [thetad], [z], [h], [theta]])
inputs = Matrix([[F], [tau]])

A = svf.jacobian(states)
B = svf.jacobian(inputs)

A_lin = simplify(A.subs([(thetad, 0.0), (zd, 0.0), (hd, 0.0), (theta, 0), (tau, 0.0), (F, g*(2*m1 + m2))]))
B_lin = simplify(B.subs([(thetad, 0.0), (zd, 0.0), (hd, 0.0), (theta, 0), (tau, 0.0), (F, g*(2*m1 + m2))]))

A_full = simplify(A.subs([(thetad, 0.0), (zd, 0.0), (hd, 0.0), (theta, 0), (tau, 0.0), (F, g*(2*m1 + m2)), (z, 5), (h, 5), (mu, 0.1), (m1, 1), (m2, 0.25), (g, 9.81), (J1, 0.0042), (d, 0.3)]))
B_full = simplify(B.subs([(thetad, 0.0), (zd, 0.0), (hd, 0.0), (theta, 0), (tau, 0.0), (F, g*(2*m1 + m2)), (z, 5), (h, 5), (mu, 0.1), (m1, 1), (m2, 0.25), (g, 9.81), (J1, 0.0042), (d, 0.3)]))

# Longitudinal
svf_Long = Matrix([[hd], [hdd_eom]])
states_Long = Matrix([[h], [hd]])
inputs_Long = Matrix([[F]])

A_Long = svf_Long.jacobian(states_Long)
B_Long = svf_Long.jacobian(inputs_Long)

A_Long = simplify(A_Long.subs([(thetad, 0.0), (zd, 0.0), (hd, 0.0), (theta, 0), (tau, 0.0), (F, g*(2*m1 + m2))]))
B_Long = simplify(B_Long.subs([(thetad, 0.0), (zd, 0.0), (hd, 0.0), (theta, 0), (tau, 0.0), (F, g*(2*m1 + m2))]))

# Lateral
svf_Lat = Matrix([[zd], [thetad], [zdd_eom], [thetadd_eom]])
states_Lat = Matrix([[z], [theta], [zd], [thetad]])
inputs_Lat = Matrix([[tau]])

A_Lat = svf_Lat.jacobian(states_Lat)
B_Lat = svf_Lat.jacobian(inputs_Lat)

A_Lat = simplify(A_Lat.subs([(thetad, 0.0), (zd, 0.0), (hd, 0.0), (theta, 0), (tau, 0.0), (F, g*(2*m1 + m2))]))
B_Lat = simplify(B_Lat.subs([(thetad, 0.0), (zd, 0.0), (hd, 0.0), (theta, 0), (tau, 0.0), (F, g*(2*m1 + m2))]))

display("Linear EOMs (A) then (B):")
display(Math(vlatex(A_full)))
display(Math(vlatex(B_full)))

display("Linear Longitudinal EOMs (A) then (B):")
display(Math(vlatex(A_Long)))
display(Math(vlatex(B_Long)))

display("Linear Lateral EOMs (A) then (B):")
display(Math(vlatex(A_Lat)))
display(Math(vlatex(B_Lat)))

# Getting T.F.
C = Matrix([[0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]])
D = Matrix([[0, 0], [0, 0], [0, 0]])                                                                                                                                                                                                        
I = eye(6)
s = symbols('s')

TF = simplify(C @ (s*I - A_lin).inv() @ B_lin + D)

display("TF:")
display(Math(vlatex(TF)))
display(Math(vlatex(simplify(TF[5]/TF[1]))))

# %% 
# Solve the block diagram
k_p, k_d = symbols('k_p, k_d')
# charEQ = 1/TF[0]
charEQ = s**2 + (k_d/(2*0.25 + 1))*s + (k_p/(2*0.25 + 1))
sVar = sp.solve(charEQ, s)
display(Math(vlatex(sVar)))

# eq1 = sVar[0]
# eq2 = sVar[1]

# eqn = Matrix([[eq1], [eq2]])

# res = sp.solve(eqn, k_d, k_p)
# display(Math(vlatex(res)))

#%%