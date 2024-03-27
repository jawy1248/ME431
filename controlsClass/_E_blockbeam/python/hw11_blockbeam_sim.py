import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics
from ctrlObserver_block import ctrlObserver
from dataPlotterObserver import dataPlotterObserver

# instantiate reference input classes
zRef = signalGenerator(0.1, 0.1, 0.15)
dRef = signalGenerator(1, 0)
nRef = signalGenerator(0, 0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()
observerPlot = dataPlotterObserver()
blockbeam = blockbeamDynamics()
ctrl = ctrlObserver()
plt.pause(5)

t = P.t_start  # time starts at t_start
y = blockbeam.h()
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        r = zRef.square(t)
        d, n = 0, 0

        u, x_hat = ctrl.update(r, y+n)
        y = blockbeam.update(u+d)
        t = t + P.Ts
    # update animation
    animation.update(blockbeam.state)
    dataPlot.update(t, r, blockbeam.state, u)
    observerPlot.update(t, blockbeam.state, x_hat, 0, 0)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
