import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics

# instantiate reference input classes
fLeft_ref = SignalGenerator(0.07, 0.5, 0.35)
fRight_ref = SignalGenerator(0.07, 0.51, 0.35)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
dyn = HummingbirdDynamics()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # set variables
        fl = fLeft_ref.sin(t)
        fr = fRight_ref.sin(t)
        force = (fl + fr)
        tau = P.d * (fl - fr)
        u = np.array([[fl], [fr]])
        dyn.update(u)
        ref = np.array([[dyn.state[0][0]], [dyn.state[1][0]], [dyn.state[2][0]]])
        t = t + P.Ts
    
    # update animation
    animation.update(t, dyn.state)
    dataPlot.update(t, dyn.state, ref, force, tau)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
