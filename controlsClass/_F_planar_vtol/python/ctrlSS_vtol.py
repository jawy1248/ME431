import numpy as np
import VTOLParam as P
import control as cnt

class ctrlSS:
    def __init__(self):
        zeta = 0.707
        A, B, C, D = P.A, P.B, P.C, P.D
        Ah, Bh, Ch, Dh = P.Ah, P.Bh, P.Ch, P.Dh

        #  tuning parameters (h)
        tr_h = 1.75        # tuned for faster rise time before saturation.
        wn_h = (0.5*np.pi) / (tr_h * np.sqrt(1 - zeta**2))
        char1 = [1, 2*zeta*wn_h, wn_h**2]
        des_poles_long = np.roots(char1)

        #  tuning parameters (theta)
        tr_theta = 0.175          # tuned for faster rise time before saturation.
        wn_theta = (0.5*np.pi) / (tr_theta * np.sqrt(1 - zeta**2))
        char2 = [1, 2*zeta*wn_theta, wn_theta**2]

        #  tuning parameters (z)
        tr_z = 10*tr_theta          # tuned for faster rise time before saturation.
        wn_z = (0.5*np.pi) / (tr_z * np.sqrt(1 - zeta**2))
        char3 = [1, 2*zeta*wn_z, wn_z**2]
        des_poles_lat = np.roots(np.convolve(char2, char3))

        # --------------- Longitudinal Dynamics ---------------
        if np.linalg.matrix_rank(cnt.ctrb(Ah, Bh)) != 2:
            print("The system is not controllable")
        else:
            self.Kh = cnt.place(Ah, Bh, des_poles_long)
            self.krh = -1.0 / (Ch @ np.linalg.inv(Ah - Bh @ self.Kh) @ Bh)
        # print gains to terminal
        print('Longitudinal Dynamics:')
        print('K: ', self.Kh)
        print('kr: ', self.krh)
        print('poles: ', des_poles_long)

        # -------------------- Lateral Dynamics --------------------
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
            print("The system is not controllable")
        else:
            self.K = cnt.place(A, B, des_poles_lat)
            Cr = np.array([[1, 0, 0, 0]])
            self.kr = -1.0 / (Cr @ np.linalg.inv(A - B @ self.K) @ B)
        # print gains to terminal
        print('Lateral Dynamics:')
        print('K: ', self.K)
        print('kr: ', self.kr)
        print('poles: ', des_poles_lat)

        self.h_d1 = 0
        self.z_d1 = 0
        self.theta_d1 = 0

        self.h_dot = 0
        self.z_dot = 0
        self.theta_dot = 0

        # dirty derivative gains
        self.sigma = 0.05
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts)

    def update(self, z_r, h_r, y):
        # x = [z, h, theta, zdot, hdot, thetadot]
        # Get states
        z = y[0][0]
        h = y[1][0]
        theta = y[2][0]
        self.z_dot = self.beta * self.z_dot + (1 - self.beta) * ((z - self.z_d1) / P.Ts)
        self.h_dot = self.beta * self.h_dot + (1 - self.beta) * ((h - self.h_d1) / P.Ts)
        self.theta_dot = self.beta * self.theta_dot + (1 - self.beta) * ((theta - self.theta_d1) / P.Ts)

        # Solve for force and torque
        stateH = np.array([[h], [self.h_dot]])
        F_unsat = -self.Kh @ stateH + self.krh * h_r
        F_unsat = F_unsat[0][0]

        state = np.array([[z], [theta], [self.z_dot], [self.theta_dot]])
        T_unsat = -self.K @ state + self.kr * z_r
        T_unsat = T_unsat[0][0]

        Fe = P.g * (2*P.mr + P.mc)
        F_unsat = F_unsat + Fe

        Fr = 0.5*F_unsat + 0.5*T_unsat/P.d
        Fl = 0.5*F_unsat - 0.5*T_unsat/P.d

        # Saturate forces
        F = np.array([[saturate(Fr, P.fmax)], [saturate(Fl, P.fmax)]])

        # update delayed variables
        self.z_d1 = z
        self.h_d1 = h
        self.theta_d1 = theta

        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








