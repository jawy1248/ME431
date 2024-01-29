import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter

# instantiate reference input classes
reference = signalGenerator(amplitude=0.5, frequency=0.2)
zRef = signalGenerator(amplitude=1.0, frequency=0.1, y_offset=2)
hRef = signalGenerator(amplitude=1.0, frequency=0.2, y_offset=2)
thetaRef = signalGenerator(amplitude=np.pi, frequency=np.pi/10)
fRef = signalGenerator(amplitude=5, frequency=.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    z = zRef.sin(t)
    h = hRef.sin(t)
    theta = thetaRef.sin(t)
    f1 = 0
    f2 = 0
    # update animation
    state = np.array([[z], [h], [theta]])  #state is made of z, h, and theta
    animation.update(state)
    dataPlot.update(t, state, z, h, f1, f2)
    # advance time by t_plot
    t = t + P.t_plot  
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
