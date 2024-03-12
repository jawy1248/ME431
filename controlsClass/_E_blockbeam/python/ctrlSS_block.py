import control as cnt
import numpy as np
import blockbeamParam as P

class ctrlSS:
    def __init__(self):
        zeta = 0.707
        A, B, C, D = P.A, P.B, P.C, P.D

        #  tuning parameters (theta)
        tr_theta = 0.1         # tuned for faster rise time before saturation.
        wn_theta = (0.5*np.pi) / (tr_theta * np.sqrt(1 - zeta**2))
        char1 = [1, 2*zeta*wn_theta, wn_theta**2]

        #  tuning parameters (z)
        tr_z = 10*tr_theta          # tuned for faster rise time before saturation.
        wn_z = (0.5*np.pi) / (tr_z * np.sqrt(1 - zeta**2))
        char2 = [1, 2*zeta*wn_z, wn_z**2]

        des_char_poly = np.convolve(char1, char2)
        des_poles = np.roots(des_char_poly)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
            print("The system is not controllable")
        else:
            self.K = cnt.place(A, B, des_poles)
            Cr = np.array([[1, 0, 0, 0]])
            self.kr = -1.0 / (Cr @ np.linalg.inv(A - B @ self.K) @ B)

        # print gains to terminal
        print('K: ', self.K)
        print('kr: ', self.kr)
        print('poles: ', des_poles)

        self.z_d1 = 0
        self.theta_d1 = 0

        self.z_dot = 0
        self.theta_dot = 0

        # dirty derivative gains
        self.sigma = 0.05
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts)

    def update(self, z_r, y):
        # Get states
        z = y[0][0]
        theta = y[1][0]
        self.z_dot = self.beta * self.z_dot + (1 - self.beta) * ((z - self.z_d1) / P.Ts)
        self.theta_dot = self.beta * self.theta_dot + (1 - self.beta) * ((theta - self.theta_d1) / P.Ts)
        state = np.array([[z], [theta], [self.z_dot], [self.theta_dot]])

        Fe = P.m2*P.g/2

        force_unsat = -self.K @ state + self.kr * z_r
        force_unsat = force_unsat[0][0] + Fe
        force = saturate(force_unsat, P.Fmax)

        # update delayed variables
        self.z_d1 = z
        self.theta_d1 = theta

        return force

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
