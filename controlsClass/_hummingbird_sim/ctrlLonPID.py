import numpy as np
import hummingbirdParam as P


class ctrlLonPID:
    def __init__(self):
        # tuning parameters
        tr_pitch = 1
        zeta_pitch = 0.707
        self.ki_pitch = 1

        # gain calculation
        b_theta = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)

        #print('b_theta: ', b_theta)
        wn_pitch = (0.5*np.pi) / (tr_pitch * np.sqrt(1 - zeta_pitch**2))

        alpha1_pitch = 2*zeta_pitch*wn_pitch
        alpha2_pitch = wn_pitch**2

        self.kd_pitch = alpha1_pitch/b_theta
        self.kp_pitch = alpha2_pitch/b_theta

        # print gains to terminal
        print('kp_pitch: ', self.kp_pitch)
        print('ki_pitch: ', self.ki_pitch)
        print('kd_pitch: ', self.kd_pitch)

        # sample rate of the controller
        self.Ts = P.Ts

        # dirty derivative parameters
        sigma = 0.01  # cutoff freq for dirty derivative
        self.beta = (2 * sigma - self.Ts) / (2 * sigma + self.Ts)

        # delayed variables
        self.theta_d1 = 0.
        self.theta_dot = 0.
        self.integrator_theta = 0.
        self.error_theta_d1 = 0.  # pitch error delayed by 1

    def update(self, r, y):
        theta_ref = r[0][0]
        theta = y[1][0]
        force_fl = P.g * (P.m1*P.ell1 + P.m2*P.ell2) / P.ellT

        # compute errors
        error_theta = theta_ref - theta

        # update differentiators
        self.theta_dot = self.beta * self.theta_dot + (1 - self.beta) * ((theta - self.theta_d1) / P.Ts)
        
        # update integrators
        self.integrator_theta = self.integrator_theta + (P.Ts / 2) * (error_theta + self.error_theta_d1)
        
        # pitch control
        force_unsat = (self.kp_pitch * error_theta) + (self.ki_pitch * self.integrator_theta) - (self.kd_pitch * self.theta_dot) + force_fl
        force = saturate(force_unsat, -P.force_max, P.force_max)
        torque = 0.

        # convert force and torque to pwm signals
        pwm = np.array([[force + torque / P.d],               # u_left
                      [force - torque / P.d]]) / (2 * P.km)   # r_right          
        pwm = saturate(pwm, 0, 1)

        # update all delayed variables
        self.theta_d1 = theta
        self.error_theta_d1 = error_theta

        # return pwm plus reference signals
        return pwm, np.array([[0.], [theta_ref], [0.]])


def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        if u > up_limit:
            u = up_limit
        if u < low_limit:
            u = low_limit
    else:
        for i in range(0, u.shape[0]):
            if u[i][0] > up_limit:
                u[i][0] = up_limit
            if u[i][0] < low_limit:
                u[i][0] = low_limit
    return u




