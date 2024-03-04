#%%
# Imports
import sympy as sp
from sympy import eye, sin, cos, diff, Matrix, symbols, simplify
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

#%%
# # Problem 3 
# m1, m2, J1, k, t, r = symbols('m1, m2, J1, k, t, r')
# theta = dynamicsymbols('theta')
# phi = dynamicsymbols('phi')
# q = Matrix([[theta], [phi]])

# P1 = Matrix([[0], [0], [0]])
# P2 = Matrix([[r * cos(theta + phi)], [r * sin(theta + phi)], [0]])

# V1 = P1.diff(t)
# V2 = P2.diff(t)

# omega = Matrix([[0], [0], [theta.diff(t)]])
# J = Matrix([[0, 0, 0], [0, 0, 0], [0, 0, J1]])

# KE = simplify(0.5*m1*V1.T*V1 + 0.5*m2*V2.T*V2 + 0.5*omega.T*J*omega)

# display(Math(vlatex(KE)))

#%%
# Problem 6
m1, m2, t, ell = symbols('m1, m2, t, ell')
theta = dynamicsymbols('theta')
x = dynamicsymbols('x')
q = Matrix([[theta], [x]])

P1 = Matrix([[x], [0], [0]])
P2 = Matrix([[x + ell*cos(theta)], [ell * sin(theta)], [0]])

V1 = P1.diff(t)
V2 = P2.diff(t)

KE = simplify(0.5*m1*V1.T*V1 + 0.5*m2*V2.T*V2)

display(Math(vlatex(simplify(KE))))

#%%