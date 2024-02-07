import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import Dynamics

# instantiate reference input classes
zRef = signalGenerator(0.5, 0.2)
hRef = signalGenerator(0.5, 0.2)
fr_Ref = signalGenerator(2.0, 0.8, 7.36)
fl_Ref = signalGenerator(2.0, 0.75, 7.34)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()
VTOL = Dynamics()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        zr = 0
        hr = 0
        f1 = fr_Ref.sin(t)
        f2 = fl_Ref.sin(t)
        force = f1 + f2
        tau = P.d*(f1 - f2)
        u = np.array([[f1],[f2]])
        VTOL.update(u)
        t = t + P.Ts
    
    # update animation
    animation.update(VTOL.state)
    dataPlot.update(t, VTOL.state, zr, hr, force, tau)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()