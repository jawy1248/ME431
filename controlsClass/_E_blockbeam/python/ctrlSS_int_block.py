import control as cnt
import numpy as np
import blockbeamParam as P

class ctrlSS:
    def __init__(self):
        zeta = 0.707
        int_pole = -5

        A, B, C, D = P.A, P.B, P.C, P.D
        Cr = np.array([[1, 0, 0, 0]])
        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A, 1), 1)))),
                        np.hstack((-Cr, np.array([[0.0]])))))
        B1 = np.vstack((B, 0))

        #  tuning parameters (theta)
        tr_theta = 0.1         # tuned for faster rise time before saturation.
        wn_theta = (0.5*np.pi) / (tr_theta * np.sqrt(1 - zeta**2))
        char1 = [1, 2*zeta*wn_theta, wn_theta**2]

        #  tuning parameters (z)
        tr_z = 10*tr_theta          # tuned for faster rise time before saturation.
        wn_z = (0.5*np.pi) / (tr_z * np.sqrt(1 - zeta**2))
        char2 = [1, 2*zeta*wn_z, wn_z**2]

        char3 = [1, -int_pole]

        char4 = np.convolve(char1, char2)
        des_char_poly = np.convolve(char3, char4)

        des_poles = np.roots(des_char_poly)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 5:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]

        # print gains to terminal
        print('K: ', self.K)
        print('ki: ', self.ki)
        print('poles: ', des_poles)

        self.integrator = 0
        self.error_d1 = 0

    def update(self, z_r, y):
        # Get states
        z = y[0][0]
        z_dot = y[1][0]
        theta = y[2][0]
        theta_dot = y[3][0]

        error_z = z_r - z
        self.integrator = self.integrator + (P.Ts / 2.0) * (error_z + self.error_d1)
        self.error_d1 = error_z

        Fe = P.m2*P.g/2

        force_unsat = -self.K @ y - self.ki * self.integrator
        force_unsat = force_unsat[0] + Fe
        force = saturate(force_unsat, P.Fmax)

        return force

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
