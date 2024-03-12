import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from ctrlLonPD import ctrlLonPD

# instantiate reference input classes
thetaRef = SignalGenerator(np.pi/32, 0.1, np.pi/8)
dRef = SignalGenerator(0.01)
nRef = SignalGenerator(0.001)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
dyn = HummingbirdDynamics()
ctrl = ctrlLonPD()
plt.pause(2)

t = P.t_start  # time starts at t_start
y = dyn.h()
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # references
        ref_theta = thetaRef.square(t)
        # d, n = dRef.random(t), nRef.random(t)
        d, n = 0, 0

        r = np.array([[ref_theta], [0]])
        pwm, ref = ctrl.update(r, y + n)

        force = (pwm[0][0] + pwm[1][0]) * P.km
        torque = (P.d * (pwm[0][0] - pwm[1][0])) * P.km

        y = dyn.update(pwm + d)
        t = t + P.Ts

    # update animation
    animation.update(t, dyn.state)
    dataPlot.update(t, dyn.state, ref, force, torque)
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
