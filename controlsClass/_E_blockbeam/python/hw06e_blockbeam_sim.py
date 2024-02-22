import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
from ctrlPD_6e import ctrlPD

# instantiate reference input classes
zRef = signalGenerator(0.15, 0.01, 0.25)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()
blockbeam = blockbeamDynamics()
ctrlPD = ctrlPD()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        zr = zRef.square(t)
        u = ctrlPD.update(zr, blockbeam.state)
        blockbeam.update(u)
        t = t + P.Ts
    # update animation
    animation.update(blockbeam.state)
    dataPlot.update(t, zr, blockbeam.state, u)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
