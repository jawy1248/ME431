import numpy as np
import massParam as P

class ctrlPD:
    def __init__(self):
        self.kp = 4.5
        self.kd = 12

    def update(self, z_r, state):
        z = state[0][0]
        zd = state[1][0]
        # compute the linearized force using PD
        Fe = P.k*z
        force = (self.kp * (z_r - z)) - (self.kd * zd) + Fe
        # compute total force
        force = saturate(force, P.F_max)
        return force


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








