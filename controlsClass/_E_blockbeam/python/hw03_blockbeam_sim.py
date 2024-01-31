import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics

# instantiate reference input classes
reference = signalGenerator(amplitude=0.5, frequency=0.2)
fRef = signalGenerator(amplitude=0.1, frequency=10, y_offset=11.496)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()
blockbeam = blockbeamDynamics()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        r = 0
        u = fRef.sin(t)
        blockbeam.update(u)
        t = t + P.Ts
    # update animation
    animation.update(blockbeam.state)
    dataPlot.update(t, r, blockbeam.state, u)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
