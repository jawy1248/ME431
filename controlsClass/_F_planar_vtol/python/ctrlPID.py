import numpy as np
import VTOLParam as P

class ctrlPID:
    def __init__(self):
        zeta = 0.707

        #  tuning parameters (h)
        tr_h = 2        # tuned for faster rise time before saturation.
        wn_h = (0.5*np.pi) / (tr_h * np.sqrt(1 - zeta**2))
        alpha1_h = 2.0 * zeta * wn_h
        alpha0_h = wn_h**2
        self.kp_h = alpha0_h*(2*P.mr + P.mc)
        self.kd_h = alpha1_h*(2*P.mr + P.mc)
        self.ki_h = 0.01
        print('kp_h: ', self.kp_h)
        print('kd_h: ', self.kd_h)
        print('ki_h: ', self.ki_h)

        #  tuning parameters (theta)
        tr_theta = 0.2          # tuned for faster rise time before saturation.
        wn_theta = (0.5*np.pi) / (tr_theta * np.sqrt(1 - zeta**2))
        alpha1_theta = 2.0 * zeta * wn_theta
        alpha0_theta = wn_theta**2
        self.kp_theta = alpha0_theta*(P.Jc + 2*P.d**2*P.mr)
        self.kd_theta = alpha1_theta*(P.Jc + 2*P.d**2*P.mr)
        print('kp_theta: ', self.kp_theta)
        print('kd_theta: ', self.kd_theta)

        #  tuning parameters (z)
        tr_z = 10*tr_theta          # tuned for faster rise time before saturation.
        wn_z = (0.5*np.pi) / (tr_z * np.sqrt(1 - zeta**2))
        alpha1_z = 2.0 * zeta * wn_z
        alpha0_z = wn_z**2
        self.kp_z = alpha0_z / (-P.g)
        self.kd_z = (alpha1_z - (P.mu/(2*P.mr + P.mc))) / (-P.g)
        self.ki_z = -0.00001
        print('kp_z: ', self.kp_z)
        print('kd_z: ', self.kd_z)
        print('ki_z: ', self.ki_z)

        self.errH_n1 = 0
        self.errZ_n1 = 0
        self.h_n1 = 0
        self.z_n1 = 0
        self.theta_n1 = 0

        self.h_dot = 0
        self.z_dot = 0
        self.theta_dot = 0

        self.intH = 0
        self.intZ = 0

        # dirty derivative gains
        self.sigma = 0.05
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts)

    def update(self, z_r, h_r, y):
        # Get states
        z = y[0][0]
        h = y[1][0]
        theta = y[2][0]

        # ------------------------- Longitudinal -------------------------
        errH = h_r - h

        # Get integrator
        if errH < 0.5:
            self.intH = self.intH + (P.Ts / 2) * (errH + self.errH_n1)
        else:
            self.intH = 0

        Fe = P.g * (2*P.mr + P.mc)
        self.h_dot = self.beta * self.h_dot + (1 - self.beta) * ((h - self.h_n1) / P.Ts)

        force = (self.kp_h * errH) + (self.ki_h * self.intH) - (self.kd_h * self.h_dot) + Fe
        # ----------------------------------------------------------------

        # ---------------------------- Lateral ----------------------------
        # ------------------------ Outer Loop ------------------------
        errZ = z_r - z

        # Get integrator
        if errZ < 0.01:
            self.intZ = self.intZ + (P.Ts / 2) * (errZ + self.errZ_n1)
        else:
            self.intZ = 0

        self.z_dot = self.beta * self.z_dot + (1 - self.beta) * ((z - self.z_n1) / P.Ts)
        theta_r = (self.kp_z * errZ) + (self.ki_z * self.intZ) - (self.kd_z * self.z_dot)
        # ------------------------------------------------------------

        # ------------------------ Inner Loop ------------------------
        errT = theta_r - theta
        self.theta_dot = self.beta * self.theta_dot + (1 - self.beta) * ((theta - self.theta_n1) / P.Ts)

        tau = (self.kp_theta * errT) - (self.kd_theta * self.theta_dot)
        # ------------------------------------------------------------
        # -----------------------------------------------------------------

        # update delayed variables
        self.errZ_n1 = errZ
        self.errH_n1 = errH
        self.z_n1 = z
        self.h_n1 = h
        self.theta_n1 = theta

        # Saturate forces
        Fl = 0.5*force - (1/(2*P.d))*tau
        Fr = 0.5*force + (1/(2*P.d))*tau
        Fl = saturate(Fl, P.fmax)
        Fr = saturate(Fr, P.fmax)
        F = [[Fr], [Fl]]

        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








