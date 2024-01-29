import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter

# instantiate reference input classes
reference = signalGenerator(amplitude=0.5, frequency=0.2)
zRef = signalGenerator(amplitude=3.0, frequency=0.2)
fRef = signalGenerator(amplitude=5, frequency=.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    r = 0
    z = zRef.sin(t)
    f = 0
    # update animation
    state = np.array([[z]])  #state is made of z
    animation.update(state)
    dataPlot.update(t, r, state, f)
    # advance time by t_plot
    t = t + P.t_plot
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
