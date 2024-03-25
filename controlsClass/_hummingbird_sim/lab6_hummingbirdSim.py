import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from ctrlPD import ctrlPD

# instantiate reference input classes
thetaRef = SignalGenerator(np.pi/6, 0.1)
psiRef = SignalGenerator(np.pi/3, 0.2)
dRef = SignalGenerator(0.01)
nRef = SignalGenerator(0.001)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
dyn = HummingbirdDynamics()
ctrl = ctrlPD()
plt.pause(4)

t = P.t_start  # time starts at t_start
y = dyn.h()
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # references
        ref_theta = thetaRef.square(t)
        ref_psi = psiRef.sin(t)
        d, n = 0, 0

        r = np.array([[ref_theta], [ref_psi]])
        pwm, ref = ctrl.update(r, y)

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
