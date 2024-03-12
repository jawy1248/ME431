import numpy as np
import massParam as P

class ctrlPID:
    def __init__(self):
        #  tuning parameters
        tr = 2  # tuned for faster rise time before saturation.
        zeta = 0.70
        # desired natural frequency
        wn = (0.5*np.pi) / (tr * np.sqrt(1 - zeta**2))
        print('wn: ', wn)
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2
        # compute PD gains
        self.kp = alpha0*P.m - P.k
        self.kd = alpha1*P.m - P.b
        self.ki = 4

        print('kp: ', self.kp)
        print('kd: ', self.kd)
        print('ki: ', self.ki)

        self.errd = 0
        self.errn1 = 0
        self.zd = 0
        self.zn1 = 0
        self.int = 0

        # dirty derivative gains
        self.sigma = 0.05
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts)

    def update(self, z_r, y):
        z = y[0][0]

        # Get linearized equilibrium force
        ze = z
        # ze = 0
        Fe = P.k * ze

        # Get error
        err = z_r - z

        # Get integrator and zdot term
        if err < 0.01:
            self.int = self.int + (P.Ts / 2) * (err + self.errn1)
        else:
            self.int = 0

        self.zd = self.beta * self.zd + (1 - self.beta) * ((z - self.zn1) / P.Ts)

        # PID Control
        force = (self.kp * err) + (self.ki * self.int) - (self.kd * self.zd) + Fe

        # compute total force
        force_sat = saturate(force, P.F_max)

        # update delayed variables
        self.errn1 = err
        self.zn1 = z

        return force_sat


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








