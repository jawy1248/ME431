import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from dataPlotterObserver import dataPlotterObserver
from VTOLDynamics import Dynamics
from ctrlObserver_vtol import ctrlObserver

# instantiate reference input classes
hRef = signalGenerator(2, 0.1, 3)
zRef = signalGenerator(2, 0.1, 3)
dRef = signalGenerator(0, 0)
nRef = signalGenerator(0, 0)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()
dataPlotObs = dataPlotterObserver()
VTOL = Dynamics()
ctrl = ctrlObserver()
plt.pause(5)

t = P.t_start  # time starts at t_start
y = VTOL.h()
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        zr, hr = zRef.square(t), hRef.square(t)
        dR, dL, n = 0, 0, 0

        u, x_hat_lat, x_hat_lon = ctrl.update(zr, hr, y)
        y = VTOL.update(u)

        F = u[0][0] + u[1][0]
        Tau = P.d * (u[0][0] - u[1][0])

        t = t + P.Ts
    
    # update animation
    animation.update(VTOL.state)
    dataPlot.update(t, VTOL.state, zr, hr, F, Tau)
    dataPlotObs.update(t, VTOL.state, x_hat_lat, x_hat_lon)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()