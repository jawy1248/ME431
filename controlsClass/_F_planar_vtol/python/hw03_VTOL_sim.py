import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import Dynamics

# instantiate reference input classes
zRef = signalGenerator(amplitude=0.5, frequency=0.2)
hRef = signalGenerator(amplitude=0.5, frequency=0.2)
fl_Ref = signalGenerator(amplitude=1, frequency=0.2, y_offset=7.3)
fr_Ref = signalGenerator(amplitude=1, frequency=0.2, y_offset=7.3)

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
        u = np.array([[f1],[f2]])
        VTOL.update(u)
        t = t + P.Ts
    
    # update animation
    animation.update(VTOL.state)
    dataPlot.update(t, VTOL.state, zr, hr, f1, f2)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()