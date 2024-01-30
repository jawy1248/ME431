import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics

# instantiate reference input classes
reference = signalGenerator(amplitude=0.5, frequency=0.2)
fRef = signalGenerator(amplitude=5, frequency=0.1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()
mass = massDynamics()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables and update dynamics
        r = 0
        u = fRef.sin(t)
        z = mass.update(u)
        t = t + P.Ts
    # update animation
    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, u)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
