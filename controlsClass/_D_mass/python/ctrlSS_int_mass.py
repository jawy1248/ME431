import numpy as np
import control as cnt
import massParam as P

class ctrlSS:
    def __init__(self):
        #  tuning parameters
        tr = 0.8  # tuned for faster rise time before saturation.
        zeta = 0.70
        int_pole = -1.4

        # desired natural frequency
        wn = (0.5*np.pi) / (tr * np.sqrt(1 - zeta**2))

        des_char_state = [1, 2*zeta*wn, wn**2]
        des_char_int = [1, -int_pole]
        des_char_poly = np.convolve(des_char_state, des_char_int)
        des_poles = np.roots(des_char_poly)

        A, B, C, D = P.A, P.B, P.C, P.D
        Cr = np.array([[1, 0]])

        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A, 1), 1)))),
                        np.hstack((-Cr, np.array([[0.0]])))))
        B1 = np.vstack( (B, 0) )

        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 3:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:2]
            self.ki = K1[0][2]
        print('K: ', self.K)
        print('ki: ', self.ki)
        print(des_poles)

        self.integrator = 0
        self.error_d1 = 0

    def update(self, z_r, y):
        z = y[0][0]
        z_dot = y[1][0]

        error_z = z_r - z
        self.integrator = self.integrator + (P.Ts / 2.0) * (error_z + self.error_d1)
        self.error_d1 = error_z

        # SS Control
        force_unsat = -self.K @ y - self.ki * self.integrator

        # compute total force
        force = saturate(force_unsat[0], P.F_max)

        return force

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


