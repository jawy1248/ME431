
import numpy as np
import blockbeamParam as P

class ctrlPID:
    def __init__(self):
        zeta = 0.707

        #  tuning parameters (theta)
        tr_theta = 0.1         # tuned for faster rise time before saturation.
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
        self.ki_z = -0.4

        print('kp_z: ', self.kp_z)
        print('kd_z: ', self.kd_z)
        print('ki_z: ', self.ki_z)

        self.err_n1 = 0
        self.z_n1 = 0
        self.theta_n1 = 0

        self.z_dot = 0
        self.theta_dot = 0

        self.int = 0

        # dirty derivative gains
        self.sigma = 0.05
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts)

    def update(self, z_r, y):
        # Get states
        z = y[0][0]
        theta = y[1][0]

        # ------------------------- Outer Loop -------------------------
        errZ = z_r - z

        # Get integrator
        if errZ < 0.01:
            self.int = self.int + (P.Ts / 2) * (errZ + self.err_n1)
        else:
            self.int = 0

        self.z_dot = self.beta * self.z_dot + (1 - self.beta) * ((z - self.z_n1) / P.Ts)
        theta_r = (self.kp_z * errZ) + (self.ki_z * self.int) - (self.kd_z * self.z_dot)
        # --------------------------------------------------------------

        # ------------------------- Inner Loop -------------------------
        errTheta = theta_r - theta
        self.theta_dot = self.beta * self.theta_dot + (1 - self.beta) * ((theta - self.theta_n1) / P.Ts)

        Fe = (P.m1 * P.g * z) / P.length + (P.m2 * P.g)/2
        force = (self.kp_theta * errTheta) - (self.kd_theta * self.theta_dot) + Fe
        force_sat = saturate(force, P.Fmax)
        # --------------------------------------------------------------

        # update delayed variables
        self.err_n1 = errZ
        self.z_n1 = z
        self.theta_n1 = theta

        return force_sat

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
