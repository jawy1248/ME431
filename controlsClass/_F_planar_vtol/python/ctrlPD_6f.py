import numpy as np
import VTOLParam as P

class ctrlPD:
    def __init__(self):
        zeta = 0.707

        #  tuning parameters (h)
        tr_h = 1        # tuned for faster rise time before saturation.
        wn_h = (0.5*np.pi) / (tr_h * np.sqrt(1 - zeta**2))
        alpha1_h = 2.0 * zeta * wn_h
        alpha0_h = wn_h**2
        self.kp_h = alpha0_h*(2*P.mr + P.mc)
        self.kd_h = alpha1_h*(2*P.mr + P.mc)
        print('kp_h: ', self.kp_h)
        print('kd_h: ', self.kd_h)

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
        print('kp_z: ', self.kp_z)
        print('kd_z: ', self.kd_z)

    def update(self, z_r, h_r, state):
        # Get states
        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        zd = state[3][0]
        hd = state[4][0]
        thetad = state[5][0]

        Fe = P.g * (2*P.mr + P.mc)
        force = (self.kp_h * (h_r - h)) - (self.kd_h * hd) + Fe
        theta_r = (self.kp_z * (z_r - z)) - (self.kd_z * zd)
        tau = (self.kp_theta * (theta_r - theta)) - (self.kd_theta * thetad)

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








