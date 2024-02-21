import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import Dynamics
from ctrlPD import ctrlPD

# instantiate reference input classes
hRef = signalGenerator(1, 0.08)
zRef = signalGenerator(2.5, 0.08, 3)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()
VTOL = Dynamics()
ctrlPD = ctrlPD()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        zr = zRef.square(t)
        hr = hRef.step(t)
        u = ctrlPD.update(zr, hr, VTOL.state)
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