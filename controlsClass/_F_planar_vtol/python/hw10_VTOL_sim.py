import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import Dynamics
from ctrlSS_int_vtol import ctrlSS

# instantiate reference input classes
hRef = signalGenerator(2, 0.1, 3)
zRef = signalGenerator(2, 0.1, 3)
dRef = signalGenerator(0, 0)
nRef = signalGenerator(0, 0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()
VTOL = Dynamics(0.2)
ctrlSS = ctrlSS()
plt.pause(5)

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        zr, hr = zRef.random(t), hRef.random(t)
        dR, dL, n = 0, 0, 0

        u = ctrlSS.update(zr, hr, VTOL.state+n)
        VTOL.update(u)

        F = u[0][0] + u[1][0]
        Tau = P.d * (u[1][0] - u[0][0])

        t = t + P.Ts
    
    # update animation
    animation.update(VTOL.state)
    dataPlot.update(t, VTOL.state, zr, hr, F, Tau)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()