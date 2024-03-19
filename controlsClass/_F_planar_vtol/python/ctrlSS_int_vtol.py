import numpy as np
import VTOLParam as P
import control as cnt

class ctrlSS:
    def __init__(self):
        zeta = 0.707
        int_pole_z = -4.0
        int_pole_h = -4.0

        A, B, C, D = P.A, P.B, P.C, P.D
        Ah, Bh, Ch, Dh = P.Ah, P.Bh, P.Ch, P.Dh

        Cr = np.array([[1, 0, 0, 0]])
        Crh = np.array([[1, 0]])

        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A, 1), 1)))),
                        np.hstack((-Cr, np.array([[0.0]])))))
        B1 = np.vstack((B, 0))

        A1h = np.vstack((np.hstack((Ah, np.zeros((np.size(Ah, 1), 1)))),
                        np.hstack((-Crh, np.array([[0.0]])))))
        B1h = np.vstack((Bh, 0))

        #  tuning parameters (h)
        tr_h = 1.75        # tuned for faster rise time before saturation.
        wn_h = (0.5*np.pi) / (tr_h * np.sqrt(1 - zeta**2))
        char1 = [1, 2*zeta*wn_h, wn_h**2]
        char2 = [1, -int_pole_h]

        des_poles_long = np.roots(np.convolve(char1, char2))

        #  tuning parameters (theta)
        tr_theta = 0.2          # tuned for faster rise time before saturation.
        wn_theta = (0.5*np.pi) / (tr_theta * np.sqrt(1 - zeta**2))
        char3 = [1, 2*zeta*wn_theta, wn_theta**2]
        char4 = [1, -int_pole_z]
        char5 = np.convolve(char3, char4)

        #  tuning parameters (z)
        tr_z = 10*tr_theta          # tuned for faster rise time before saturation.
        wn_z = (0.5*np.pi) / (tr_z * np.sqrt(1 - zeta**2))
        char6 = [1, 2*zeta*wn_z, wn_z**2]
        des_poles_lat = np.roots(np.convolve(char5, char6))

        # --------------- Longitudinal Dynamics ---------------
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

        # -------------------- Lateral Dynamics --------------------
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

        self.integrator_z = 0
        self.integrator_h = 0

        self.error_d1_z = 0
        self.error_d1_h = 0

    def update(self, z_r, h_r, y):
        # x = [z, h, theta, zdot, hdot, thetadot]
        # Get states
        z = y[0][0]
        h = y[1][0]
        theta = y[2][0]
        z_dot = y[3][0]
        h_dot = y[4][0]
        theta_dot = y[5][0]

        # Solve for force and torque
        error_h = h_r - h
        self.integrator_h = self.integrator_h + (P.Ts / 2.0) * (error_h + self.error_d1_h)
        self.error_d1_h = error_h
        Fe = P.g * (2 * P.mr + P.mc)
        long_y = np.array([[h], [h_dot]])
        F_unsat = -self.Kh @ long_y - self.kih * self.integrator_h
        F_unsat = F_unsat[0] + Fe

        error_z = z_r - z
        self.integrator_z = self.integrator_z + (P.Ts / 2.0) * (error_z + self.error_d1_z)
        self.error_d1_z = error_z
        lat_y = np.array([[z], [theta], [z_dot], [theta_dot]])
        T_unsat = -self.K @ lat_y - self.ki * self.integrator_z
        T_unsat = T_unsat[0]

        Fr = 0.5*F_unsat + 0.5*T_unsat/P.d
        Fl = 0.5*F_unsat - 0.5*T_unsat/P.d

        # Saturate forces
        F = np.array([[saturate(Fr, P.fmax)], [saturate(Fl, P.fmax)]])

        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u








