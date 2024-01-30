import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter

# instantiate reference input classes
phi_ref = SignalGenerator(amplitude=1.5, frequency=0.2)
theta_ref = SignalGenerator(amplitude=0.5, frequency=0.3)
psi_ref = SignalGenerator(amplitude=0.5, frequency=0.1)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    '''
    - roll (phi): '+' is angled to the right (from perspective of stand) and '-' is angled left
        if the eq is increasing it is moving cw (from perspective of stand) and ccw as it decreases
    - pitch (theta): '+' is up, '-' is down
        if the eq is increasing it is moving up and down as it decreases
    - yaw (psi): '+' is angled cw (from perspective of top), '-' is angled ccw (from the positive x line)
        if the eq is increasing it is moving cw (from perspective of top) and ccw as it decreases
    '''
    phi = phi_ref.sin(t)
    theta = theta_ref.sin(t)
    psi = psi_ref.sin(t)
    # update animation
    state = np.array([[phi], [theta], [psi], [0.0], [0.0], [0.0]])
    ref = np.array([[0], [0], [0]])
    force = 0
    torque = 0
    animation.update(t, state)
    dataPlot.update(t, state, ref, force, torque)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
