import control as cnt
import numpy as np
import blockbeamParam as P

class ctrlObserver:
    def __init__(self):
        # ******************** Tunable Params ********************
        zeta = 0.707
        int_pole = -5

        tr_theta = 0.12

        M_z = 10
        MO_z = 10
        MO_t = 10
        # ********************************************************

        tr_z = M_z*tr_theta
        trO_theta = tr_theta/MO_t
        trO_z = tr_z/MO_z

        self.A, self.B, self.C, self.D = P.A, P.B, P.C, P.D
        Cr = np.array([[1, 0, 0, 0]])
        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A, 1), 1)))),
                        np.hstack((-Cr, np.array([[0.0]])))))
        B1 = np.vstack((self.B, 0))

        #  tuning parameters (theta)
        wn_theta = (0.5*np.pi) / (tr_theta * np.sqrt(1 - zeta**2))
        char1 = [1, 2*zeta*wn_theta, wn_theta**2]

        #  tuning parameters (z)
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

        # Observer Design
        wnO_z = (0.5*np.pi) / (trO_z * np.sqrt(1 - zeta**2))
        wnO_theta = (0.5*np.pi) / (trO_theta * np.sqrt(1 - zeta**2))

        desO_charZ = [1, 2 * zeta * wnO_z, wnO_z ** 2]
        desO_charT = [1, 2 * zeta * wnO_theta, wnO_theta ** 2]

        desO_char = np.convolve(desO_charZ, desO_charT)
        desO_poles = np.roots(desO_char)

        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 4:
            print("The system is not observerable")
        else:
            self.L = cnt.place(self.A.T, self.C.T, desO_poles).T
        print('L^T: ', self.L.T)

        self.integrator = 0.0
        self.error_d1 = 0.0
        self.x_hat = np.array([[0.0], [0.0], [0.0], [0.0]])
        self.force_d1 = 0.0

    def update(self, z_r, y):
        # Get states
        x_hat = self.update_observer(y)
        z_hat = x_hat[0][0]

        error_z = z_r - z_hat
        self.integrator = self.integrator + (P.Ts / 2.0) * (error_z + self.error_d1)
        self.error_d1 = error_z

        Fe = P.m2 * P.g / 2
        force_unsat = -self.K @ x_hat - self.ki * self.integrator
        force = saturate(Fe + force_unsat[0], P.Fmax)
        self.force_d1 = force

        return force, x_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.x_hat, y_m)
        F2 = self.observer_f(self.x_hat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.x_hat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.x_hat + P.Ts * F3, y_m)
        self.x_hat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat

    def observer_f(self, x_hat, y_m):
        z_hat = x_hat[0][0]
        Fe = P.m2*P.g / 2
        xhat_dot = self.A @ x_hat + self.B * (self.force_d1 - Fe) + self.L @ (y_m - self.C @ x_hat)
        return xhat_dot

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
