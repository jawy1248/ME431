import numpy as np
import massParam as P

class ctrlPD:
    def __init__(self):
        #  tuning parameters
        tr = 1.5  # tuned for faster rise time before saturation.
        zeta = 0.70
        # desired natural frequency
        wn = (0.5*np.pi) / (tr * np.sqrt(1 - zeta**2))
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2
        # compute PD gains
        self.kp = alpha0*P.m - P.k
        self.kd = alpha1*P.m - P.b
        print('kp: ', self.kp)
        print('kd: ', self.kd)    

    def update(self, z_r, state):
        z = state[0][0]
        zd = state[1][0]
        # compute the linearized force using PD
        force = (self.kp * (z_r - z)) - (self.kd * zd)
        # compute total force
        force = saturate(force, P.F_max)
        return force


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








