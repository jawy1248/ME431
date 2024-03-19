import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
from ctrlSS_int_block import ctrlSS

# instantiate reference input classes
zRef = signalGenerator(0.1, 0.1, 0.15)
dRef = signalGenerator(1, 0)
nRef = signalGenerator(0, 0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()
blockbeam = blockbeamDynamics(0.2)
ctrlSS = ctrlSS()
plt.pause(5)

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        r = zRef.square(t)
        d, n = dRef.step(t), 0
        # d, n = 0, 0

        u = ctrlSS.update(r, blockbeam.state+n)
        y = blockbeam.update(u+d)
        t = t + P.Ts
    # update animation
    animation.update(blockbeam.state)
    dataPlot.update(t, r, blockbeam.state, u)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
