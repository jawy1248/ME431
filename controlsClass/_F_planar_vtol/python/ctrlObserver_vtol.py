import numpy as np
import VTOLParam as P
import control as cnt

class ctrlObserver:
    def __init__(self):
        # ******************** Tunable Params ********************
        zeta = 0.707
        int_pole_z = -3.0
        int_pole_h = -3.0

        tr_h = 2
        tr_t = 0.2

        M_z = 10
        MO_z = 30
        MO_h = 20
        MO_t = 20
        # ********************************************************

        tr_z = M_z*tr_t
        trO_z = MO_z*tr_z
        trO_h = MO_h*tr_h
        trO_t = MO_t*tr_t

        # --------------- Longitudinal Dynamics ---------------
        self.Ah, self.Bh, self.Ch, self.Dh = P.Ah, P.Bh, P.Ch, P.Dh
        Crh = np.array([[1, 0]])

        A1h = np.vstack((np.hstack((self.Ah, np.zeros((np.size(self.Ah, 1), 1)))),
                        np.hstack((-Crh, np.array([[0.0]])))))
        B1h = np.vstack((self.Bh, 0))

        wn_h = (0.5*np.pi) / (tr_h * np.sqrt(1 - zeta**2))
        char1 = [1, 2*zeta*wn_h, wn_h**2]
        char2 = [1, -int_pole_h]
        des_poles_long = np.roots(np.convolve(char1, char2))

        if np.linalg.matrix_rank(cnt.ctrb(A1h, B1h)) != 3:
            print("The system is not controllable")
        else:
            K1h = cnt.place(A1h, B1h, des_poles_long)
            self.Kh = K1h[0][0:2]
            self.kih = K1h[0][2]
        # print gains to terminal
        print('Longitudinal Dynamics:')
        print('K: ', self.Kh)
        print('ki: ', self.kih)
        print('poles: ', des_poles_long)

        # Observer Design
        wnO_h = (0.5*np.pi) / (trO_h * np.sqrt(1 - zeta**2))
        desO_charH = [1, 2 * zeta * wnO_h, wnO_h**2]
        desO_polesH = np.roots(desO_charH)

        if np.linalg.matrix_rank(cnt.ctrb(self.Ah.T, self.Ch.T)) != 2:
            print("The Longitudinal system is not observerable")
        else:
            self.Lh = cnt.place(self.Ah.T, self.Ch.T, desO_polesH).T
        print('Lh^T: ', self.Lh.T)

        # -------------------- Lateral Dynamics --------------------
        self.A, self.B, self.C, self.D = P.A, P.B, P.C, P.D
        Cr = np.array([[1, 0, 0, 0]])

        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A, 1), 1)))),
                        np.hstack((-Cr, np.array([[0.0]])))))
        B1 = np.vstack((self.B, 0))

        wn_t = (0.5*np.pi) / (tr_t * np.sqrt(1 - zeta**2))
        char3 = [1, 2*zeta*wn_t, wn_t**2]
        char4 = [1, -int_pole_z]
        char5 = np.convolve(char3, char4)

        wn_z = (0.5*np.pi) / (tr_z * np.sqrt(1 - zeta**2))
        char6 = [1, 2*zeta*wn_z, wn_z**2]
        des_poles_lat = np.roots(np.convolve(char5, char6))

        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 5:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles_lat)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]
        # print gains to terminal
        print('Lateral Dynamics:')
        print('K: ', self.K)
        print('ki: ', self.ki)
        print('poles: ', des_poles_lat)

        # Observer Design
        wnO_z = (0.5*np.pi) / (trO_z * np.sqrt(1 - zeta**2))
        wnO_t = (0.5*np.pi) / (trO_t * np.sqrt(1 - zeta**2))

        desO_charZ = [1, 2 * zeta * wnO_z, wnO_z ** 2]
        desO_charT = [1, 2 * zeta * wnO_t, wnO_t ** 2]

        desO_char = np.convolve(desO_charZ, desO_charT)
        desO_poles = np.roots(desO_char)

        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 4:
            print("The Lateral system is not observerable")
        else:
            self.L = cnt.place(self.A.T, self.C.T, desO_poles).T
        print('L^T: ', self.L.T)

        # -------------------- Misc Variables --------------------
        self.integrator_z = 0
        self.integrator_h = 0

        self.error_d1_z = 0
        self.error_d1_h = 0

        self.x_hat_long = np.array([[0.0], [0.0]])
        self.x_hat_lat = np.array([[0.0], [0.0], [0.0], [0.0]])

        self.force_d1 = 0.0
        self.torque_d1 = 0.0

    def update(self, z_r, h_r, y):
        # x_hat = [z, h, theta]
        # Get states
        x_hat_long = self.update_observer_long(y[1][0])
        x_hat_lat = self.update_observer_lat(y[0:3:2][0])

        h_hat = x_hat_long[0][0]
        z_hat = x_hat_lat[0][0]

        # Solve for force
        error_h = h_r - h_hat
        self.integrator_h = self.integrator_h + (P.Ts / 2.0) * (error_h + self.error_d1_h)
        self.error_d1_h = error_h

        Fe = P.g * (2 * P.mr + P.mc)
        F_unsat = -self.Kh @ x_hat_long - self.kih * self.integrator_h
        F_unsat = F_unsat[0] + Fe

        # Solve for torque
        error_z = z_r - z_hat
        self.integrator_z = self.integrator_z + (P.Ts / 2.0) * (error_z + self.error_d1_z)
        self.error_d1_z = error_z

        T_unsat = -self.K @ x_hat_lat - self.ki * self.integrator_z
        T_unsat = T_unsat[0]
        # T_unsat = 0

        # Combine forces and torques, saturate and return
        Fr_unsat = 0.5*F_unsat + 0.5*T_unsat/P.d
        Fl_unsat = 0.5*F_unsat - 0.5*T_unsat/P.d

        Fr = saturate(Fr_unsat, P.fmax)
        Fl = saturate(Fl_unsat, P.fmax)

        F = np.array([[Fr], [Fl]])

        self.force_d1 = Fr + Fl
        self.torque_d1 = (Fr - Fl) * P.d

        return F, x_hat_lat, x_hat_long

    def update_observer_long(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_long(self.x_hat_long, y_m)
        F2 = self.observer_f_long(self.x_hat_long + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f_long(self.x_hat_long + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f_long(self.x_hat_long + P.Ts * F3, y_m)
        self.x_hat_long += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat_long

    def observer_f_long(self, x_hat, y_m):
        Fe = P.g * (2 * P.mr + P.mc)
        xhat_dot = self.Ah @ x_hat + self.Bh * (self.force_d1 - Fe) + self.Lh * (y_m - self.Ch @ x_hat)
        return xhat_dot

    def update_observer_lat(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lat(self.x_hat_lat, y_m)
        F2 = self.observer_f_lat(self.x_hat_lat + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lat(self.x_hat_lat + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lat(self.x_hat_lat + P.Ts * F3, y_m)
        self.x_hat_lat += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        return self.x_hat_lat

    def observer_f_lat(self, x_hat, y_m):
        xhat_dot = self.A @ x_hat + self.B * self.torque_d1 + self.L @ (y_m - self.C @ x_hat)
        return xhat_dot

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








