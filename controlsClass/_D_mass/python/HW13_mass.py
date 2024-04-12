# ************************ HW13 ************************
from control import bode, tf, margin
import matplotlib.pyplot as plt
import numpy as np
import loopshape_tools as lt
import massParam as P
from discreteFilter import discreteFilter

# ************************ D.17 ************************
# Define Plant
Plant = tf([1/P.m],[1, P.b/P.m, P.k/P.m])

# Define PID Controller
kd = 7.277
kp = 3.05
ki = 1.5
C = tf([kd, kp, ki], [1.0, 0])

'''
# ************************ D.18 ************************
#integral control
k_I = 0.4
Integrator = lt.get_control_integral(k_I)
C = C*Integrator

#proportional control
kp = 1/0.77  # 1/(current Mag. ratio at desired crossover freq.)
C = C*kp

#phase lag
z = 0.8
M = 100.0
# Lag = tf([1, z], [1, z/M])
Lag = lt.get_control_lag(z, M)
C = C*Lag

#low-pass filter
p = 5
#LPF = tf([p], [1, p])
LPF = lt.get_control_lpf(p)
C = C*LPF


# this shows how to now turn our continuous-time controller transfer function
# into a discrete-time controller/filter. This is for a sampling time of 0.01 sec. 
controller = discreteFilter(C.num, C.den, 0.01)

# to use this function, it needs to be included or defined in a controller and will 
# return the desired output when called as follows:
#   u = controller.update(error)
'''

if __name__=="__main__":
    dB_flag = False

    # -------------------- Plotting --------------------
    # Plant
    mag, phase, omega = bode(Plant,
                             dB=dB_flag,
                             omega=np.logspace(-4, 5),
                             label='P(s) - Plant',
                             plot=True)

    gm, pm, Wcg, Wcp = margin(Plant)
    print("for original system:")
    print(" pm: ", pm, " Wcp: ", Wcp, "gm: ", gm, " Wcg: ", Wcg)

    # -------------------- Design Specifications --------------------
    # ----------- Reference Tracking -----------
    omega_r = 0.1
    gamma_r = 0.03
    lt.add_spec_ref_tracking(gamma_r, omega_r)

    # ----------- Noise Attenuation -----------
    omega_n = 500
    gamma_n = 0.001
    lt.add_spec_noise(gamma_n, omega_n)

    # ----------- Extra Requirements -----------
    '''
    - Reject constant disturbances (make it system type = 1)
    - Make PM close to 60 degrees
    - Use prefilter to reduce any peaking in closed-loop solution
    '''

    # -------------------- Plotting --------------------
    # Plant + Controller
    mag, phase, omega = bode(Plant*C,
                             dB=dB_flag,
                             omega=np.logspace(-4, 5),
                             label='P(s)C(s) - Open Loop',
                             plot=True,
                             margins=True)

    gm, pm, Wcg, Wcp = margin(Plant * C)
    print("for final C*P:")
    print(" pm: ", pm, " Wcp: ", Wcp, "gm: ", gm, " Wcg: ", Wcg)

    # # Closed-Loop Controller
    # mag, phase, omega = bode((Plant*C / (1 + Plant*C)),
    #                          dB=dB_flag,
    #                          omega=np.logspace(-4, 5),
    #                          label=r'$\frac{P(s)C(s)}{1+P(s)C(s)}$ - Closed-loop',
    #                          plot=True,
    #                          margins=True)

    fig = plt.gcf()
    fig.axes[0].legend()
    fig.axes[0].grid(True)
    fig.axes[1].grid(True)
    fig.axes[0].set_title('Bode Diagram')
    # plt.savefig('D18.png')
    plt.show()
