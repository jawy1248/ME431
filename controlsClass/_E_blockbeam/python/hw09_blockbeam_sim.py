import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
from ctrlSS_block import ctrlSS

# instantiate reference input classes
zRef = signalGenerator(0.05, 1)
dRef = signalGenerator(0.003, 0.08, 0.006)
nRef = signalGenerator(0.005, 0.1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()
blockbeam = blockbeamDynamics()
ctrlSS = ctrlSS()
plt.pause(2)

t = P.t_start  # time starts at t_start
y = blockbeam.h()
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        r = zRef.step(t)
        # d, n = dRef.random(t), nRef.random(t)
        d, n = 0, 0

        u = ctrlSS.update(r, y+n)
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
