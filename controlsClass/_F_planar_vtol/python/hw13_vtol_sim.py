import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import Dynamics
from ctrlLoopshapeLong import ctrlLoopshapeLong
from ctrlLoopshapeLat import ctrlLoopshapeLat

# instantiate reference input classes
zRef = signalGenerator(1, 0.1/(2*np.pi))
hRef = signalGenerator(1, 0.1/(2*np.pi))

RdRef = signalGenerator(0.25, 0)
LdRef = signalGenerator(0.25, 0)

ZnRef = signalGenerator(0.001, 500/(2*np.pi))
HnRef = signalGenerator(0.001, 500/(2*np.pi))

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()
vtol = Dynamics()
ctrlLong = ctrlLoopshapeLong(method="digital_filter")
ctrlLat = ctrlLoopshapeLat(method="digital_filter")
plt.pause(2)

t = P.t_start  # time starts at t_start
y = vtol.h()
while t < P.t_end:  # main simulation loop
    # Progress dynamics between plot samples
    t_next_plot = t + P.t_plot
    # update controls and dynamics
    while t < t_next_plot:
        # Ref
        h_ref = hRef.square(t)
        z_ref = zRef.square(t)

        # Disturbance
        rD = RdRef.step(t)
        lD = LdRef.step(t)
        d = np.array([[rD], [lD]])

        # Noise
        zN = ZnRef.random(t)
        hN = ZnRef.random(t)
        n = np.array([[zN], [hN], [0]])

        Force = ctrlLong.update(h_ref, y+n)
        Torque = ctrlLat.update(z_ref, y+n)

        fr = 0.5*Force + Torque/(2*P.d)
        fl = 0.5*Force - Torque/(2*P.d)

        u = np.array([[fr], [fl]])

        y = vtol.update(u+d)
        t = t + P.Ts

    # update animation
    animation.update(vtol.state)
    dataPlot.update(t, vtol.state, z_ref, h_ref, Force, Torque)
    plt.pause(0.001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
