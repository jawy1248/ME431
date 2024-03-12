import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import Dynamics
from ctrlSS_vtol import ctrlSS

# instantiate reference input classes
hRef = signalGenerator(2, 0.08, 4)
zRef = signalGenerator(2.5, 0.08, 3)
dRef = signalGenerator(0.3, 0.08, 0.1)
nRef = signalGenerator(0.05, 0.1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()
VTOL = Dynamics()
ctrlSS = ctrlSS()
plt.pause(2)

t = P.t_start  # time starts at t_start
y = VTOL.h()
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        zr, hr = zRef.step(t), hRef.step(t)
        # dR, dL, n = dRef.random(t), dRef.random(t), nRef.random(t)
        dR, dL, n = 0, 0, 0

        u = ctrlSS.update(zr, hr, y+n)
        if dL != 0 and dR != 0:
            u[0][0] = u[0][0] + dR
            u[1][0] = u[1][0] + dL
        y = VTOL.update(u)

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