import numpy as np
import control as cnt
import massParam as P

class ctrlSS:
    def __init__(self):
        #  tuning parameters
        tr = 2  # tuned for faster rise time before saturation.
        zeta = 0.70
        # desired natural frequency
        wn = (0.5*np.pi) / (tr * np.sqrt(1 - zeta**2))

        des_char_poly = [1, 2*zeta*wn, wn**2]
        des_poles = np.roots(des_char_poly)

        A, B, C, D = P.A, P.B, P.C, P.D
        
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 2:
            print("The system is not controllable")
        else:
            self.K = (cnt.place(A, B, des_poles))
            self.kr = -1.0 / (C @ np.linalg.inv(A - B @ self.K) @ B)
        print('K: ', self.K)
        print('kr: ', self.kr)
        print(des_poles)

        # discrete time variables
        self.zdot = 0
        self.z_d1 = 0

        # dirty derivative gains
        self.sigma = 0.05
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts)

    def update(self, z_r, y):
        z = y[0][0]
        self.zdot = self.beta * self.zdot + (1 - self.beta) * ((z - self.z_d1) / P.Ts)

        # SS Control
        state = np.array([[z], [self.zdot]])
        force_unsat = -self.K @ state + self.kr * z_r

        # compute total force
        force = saturate(force_unsat[0][0], P.F_max)

        # update delayed variables
        self.z_d1 = z

        return force

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








