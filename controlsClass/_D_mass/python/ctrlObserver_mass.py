import numpy as np
import control as cnt
import massParam as P

class ctrlObserver:
    def __init__(self):
        # ******************** Tunable Params ********************
        zeta = 0.707
        int_pole = -1.4
        tr = 0.8
        # ********************************************************

        tr_obs = tr/10

        # desired natural frequency
        wn = (0.5*np.pi) / (tr * np.sqrt(1 - zeta**2))

        des_char_state = [1, 2*zeta*wn, wn**2]
        des_char_int = [1, -int_pole]
        des_char_poly = np.convolve(des_char_state, des_char_int)
        des_poles = np.roots(des_char_poly)

        self.A, self.B, self.C, self.D = P.A, P.B, P.C, P.D
        Cr = np.array([[1, 0]])

        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A, 1), 1)))),
                        np.hstack((-Cr, np.array([[0.0]])))))
        B1 = np.vstack( (self.B, 0) )

        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 3:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:2]
            self.ki = K1[0][2]
        print('K: ', self.K)
        print('ki: ', self.ki)
        print(des_poles)

        # Observer Design
        wn_obs = (0.5*np.pi) / (tr_obs * np.sqrt(1 - zeta**2))
        des_obs_char = [1, 2*zeta*wn_obs, wn_obs**2]
        des_obs_poles = np.roots(des_obs_char)
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 2:
            print("The system is not observerable")
        else:
            self.L = cnt.place(self.A.T, self.C.T, des_obs_poles).T
        print('L^T: ', self.L.T)

        self.integrator = 0.0
        self.error_d1 = 0.0
        self.x_hat = np.array([[0.0], [0.0]])
        self.force_d1 = 0.0

    def update(self, z_r, y):
        x_hat = self.update_observer(y)
        z_hat = x_hat[0][0]

        error_z = z_r - z_hat
        self.integrator = self.integrator + (P.Ts / 2.0) * (error_z + self.error_d1)
        self.error_d1 = error_z

        # compute total force
        force_unsat = -self.K @ x_hat - self.ki * self.integrator
        force = saturate(force_unsat[0], P.F_max)
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
        z = y_m[0][0]
        xhat_dot = self.A @ x_hat + self.B * self.force_d1 + self.L * (y_m - self.C @ x_hat)
        return xhat_dot

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


