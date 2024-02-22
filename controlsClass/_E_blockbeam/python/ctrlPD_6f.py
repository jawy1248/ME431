
import numpy as np
import blockbeamParam as P

class ctrlPD:
    def __init__(self):
        zeta = 0.707

        #  tuning parameters (theta)
        tr_theta = 0.11         # tuned for faster rise time before saturation.
        wn_theta = (0.5*np.pi) / (tr_theta * np.sqrt(1 - zeta**2))
        alpha1_theta = 2.0 * zeta * wn_theta
        alpha0_theta = wn_theta**2
        self.kp_theta = alpha0_theta*((1/4)*P.m1*P.length + (1/3)*P.m2*P.length)
        self.kd_theta = alpha1_theta*((1/4)*P.m1*P.length + (1/3)*P.m2*P.length)
        print('kp_theta: ', self.kp_theta)
        print('kd_theta: ', self.kd_theta)

        #  tuning parameters (z)
        tr_z = 10*tr_theta          # tuned for faster rise time before saturation.
        wn_z = (0.5*np.pi) / (tr_z * np.sqrt(1 - zeta**2))
        alpha1_z = 2.0 * zeta * wn_z
        alpha0_z = wn_z**2
        self.kp_z = alpha0_z / (-P.g)
        self.kd_z = alpha1_z / (-P.g) 
        print('kp_z: ', self.kp_z)
        print('kd_z: ', self.kd_z)

    def update(self, z_r, state):
        # Get states
        z = state[0][0]
        theta = state[1][0]
        zd = state[2][0]
        thetad = state[3][0]

        Fe = (P.m1*P.g*z)/P.length + (P.m2*P.g)/2

        theta_r = (self.kp_z * (z_r - z)) - (self.kd_z * zd)
        force = (self.kp_theta * (theta_r - theta)) - (self.kd_theta * thetad) + Fe
        force = saturate(force, P.Fmax)

        return force


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
