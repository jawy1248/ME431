import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlLoopshape import ctrlLoopshape
from dataPlotterObserver import dataPlotterObserver

# instantiate reference input classes
zRef = signalGenerator(1, 0.1/(2*np.pi))
dRef = signalGenerator(0.25, 0)
nRef = signalGenerator(0.001, 500/(2*np.pi))

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()
mass = massDynamics()
ctrl = ctrlLoopshape(method="digital_filter")
plt.pause(2)

t = P.t_start  # time starts at t_start
y = mass.h()
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables and update dynamics
        r = zRef.step(t)
        d, n = dRef.step(t), nRef.random(t)

        u = ctrl.update(r, y+n)
        y = mass.update(u+d)
        t = t + P.Ts

    # update animation
    animation.update(mass.state)
    dataPlot.update(t, r, mass.state, u)
    plt.pause(0.001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
