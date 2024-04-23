# ************************ HW13 ************************
from control import tf, bode, margin, step_response, mag2db
import matplotlib.pyplot as plt
import numpy as np
import loopshape_tools as lt
import massParam as P

# ************************ D.17 ************************
# Define Plant
Plant = tf([1/P.m],[1, P.b/P.m, P.k/P.m])

# Define PID Controller
kd = 7.277
kp = 3.05
ki = 1.5
C = tf([kd, kp, ki], [1.0, 0])

# ************************ D.18 ************************
'''
1) get_control_integral(ki):            controls system type at low-freq
2) get_control_proportional(kp):        controls w_co
3) get_control_lead(omega_lead, M):     controls PM
4) get_control_lpf(p):                  controls noise reduction - also does prefilter
5) get_control_lag(z, M):               controls ref tracking
'''

# # # Integral control
# # Ki = 1
# # C_int = lt.get_control_integral(Ki)
# # C = C*C_int

# # Proportional control
# Kp = 10
# C_prop = lt.get_control_proportional(Kp)
# C = C*C_prop

# # Lead Compensator
# omega_z = 7
# M = 10
# C_lead = lt.get_control_lead(omega_z, M)
# C = C*C_lead

# # Low-pass filter
# omega_lpf = 1
# C_lpf = lt.get_control_lpf(omega_lpf)
# C = C*C_lpf

# # # Phase Lag
# # omega_z = 10
# # M = 10
# # C_lag = lt.get_control_lag(omega_z, M)
# # C = C*C_lag

# Get Prefilter
omega_pre = 1
F = lt.get_control_lpf(omega_pre)

# Get TF num and den
C_num = np.asarray(C.num[0])
C_den = np.asarray(C.den[0])
F_num = np.asarray(F.num[0])
F_den = np.asarray(F.den[0])


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
    # ---------------------------------------------------------------

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

    # # Closed-Loop Controller (D17)
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

    # ############################################
    # # now check the closed-loop response with prefilter
    # ############################################
    # # Closed loop transfer function from R to Y - no prefilter
    # CLOSED_R_to_Y = (Plant * C / (1.0 + Plant * C))
    # # Closed loop transfer function from R to Y - with prefilter
    # CLOSED_R_to_Y_with_F = (F * Plant * C / (1.0 + Plant * C))
    # # Closed loop transfer function from R to U - no prefilter
    # CLOSED_R_to_U = (C / (1.0 + Plant * C))
    # # Closed loop transfer function from R to U - with prefilter
    # CLOSED_R_to_U_with_F = (F*C / (1.0 + Plant * C))
    #
    # plt.figure(4)
    # plt.clf()
    # plt.grid(True)
    # plt.subplot(311)
    # mag, phase, omega = bode(CLOSED_R_to_Y, dB=dB_flag, plot=False)
    # if dB_flag:
    #     plt.semilogx(omega, mag2db(mag), color=[0, 0, 1],
    #         label='closed-loop $\\frac{Y}{R}$ - no pre-filter')
    # else:
    #     plt.loglog(omega, mag, color=[0, 0, 1],
    #         label='closed-loop $\\frac{Y}{R}$ - no pre-filter')
    # mag, phase, omega = bode(CLOSED_R_to_Y_with_F,
    #                          dB=dB_flag, plot=False)
    # if dB_flag:
    #     plt.semilogx(omega, mag2db(mag), color=[0, 1, 0],
    #         label='closed-loop $\\frac{Y}{R}$ - with pre-filter')
    # else:
    #     plt.loglog(omega, mag, color=[0, 1, 0],
    #         label='closed-loop $\\frac{Y}{R}$ - with pre-filter')
    # plt.ylabel('Closed-Loop Bode Plot')
    # plt.grid(True)
    # plt.legend()
    #
    # plt.subplot(312), plt.grid(True)
    # T = np.linspace(0, 2, 100)
    # _, yout_no_F = step_response(CLOSED_R_to_Y, T)
    # _, yout_F = step_response(CLOSED_R_to_Y_with_F, T)
    # plt.plot(T, yout_no_F, color=[0,0,1],
    #          label='response without prefilter')
    # plt.plot(T, yout_F, color=[0,1,0],
    #          label='response with prefilter')
    # plt.legend()
    # plt.ylabel('Step Response')
    #
    #
    # plt.subplot(313)
    # plt.grid(True)
    # _, Uout = step_response(CLOSED_R_to_U, T)
    # _, Uout_F = step_response(CLOSED_R_to_U_with_F, T)
    # plt.plot(T, Uout, color=[0, 0, 1],
    #          label='control effort without prefilter')
    # plt.plot(T, Uout_F, color=[0, 1, 0],
    #          label='control effort with prefilter')
    # plt.ylabel('Control Effort')
    # plt.legend()
    #
    # plt.show()
