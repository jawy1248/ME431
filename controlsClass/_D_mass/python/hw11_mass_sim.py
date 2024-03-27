import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlObserver_mass import ctrlObserver
from dataPlotterObserver import dataPlotterObserver

# instantiate reference input classes
zRef = signalGenerator(1, 0.1)
dRef = signalGenerator(0.25, 0)
nRef = signalGenerator(0, 0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()
observerPlot = dataPlotterObserver()
mass = massDynamics()
ctrl = ctrlObserver()
plt.pause(5)

t = P.t_start  # time starts at t_start
y = mass.h()
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables and update dynamics
        r = zRef.square(t)
        d, n = 0, 0

        u, x_hat = ctrl.update(r, y+n)
        y = mass.update(u+d)
        t = t + P.Ts
    # update animation
    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, u)
    observerPlot.update(t, mass.state, x_hat)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
