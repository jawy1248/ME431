import numpy as np
import hummingbirdParam as P

"""
Phi - Roll
Theta - Pitch
Psi - Yaw
"""

class ctrlPID:
    def __init__(self):
        # ******************** Tuning Parameters ********************
        # ----- Pitch -----
        tr_pitch = 1
        zeta_pitch = 0.707
        self.ki_pitch = 1

        # ----- Roll/Yaw -----
        M = 10
        tr_yaw = 1
        tr_roll = M * tr_yaw
        zeta_yaw = 0.707
        zeta_roll = 0.707
        self.ki_yaw = 1
        # ***********************************************************

        # ********************** Finding Gains **********************
        # ----- Pitch -----
        b_theta = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)

        wn_pitch = (0.5*np.pi) / (tr_pitch * np.sqrt(1 - zeta_pitch**2))

        alpha1_pitch = 2 * zeta_pitch * wn_pitch
        alpha2_pitch = wn_pitch**2

        self.kd_pitch = alpha1_pitch / b_theta
        self.kp_pitch = alpha2_pitch / b_theta

        print('kd_pitch: ', self.kd_pitch)
        print('kp_pitch: ', self.kp_pitch)
        print('ki_pitch: ', self.ki_pitch)

        # ----- Roll/Yaw -----
        self.Fe = P.g * (P.m1 * P.ell1 + P.m2 * P.ell2) / P.ellT
        b_yaw = (P.ellT * self.Fe) / (P.JT + P.J1z)

        wn_roll = (0.5*np.pi) / (tr_roll * np.sqrt(1 - zeta_roll**2))
        wn_yaw = (0.5*np.pi) / (tr_yaw * np.sqrt(1 - zeta_yaw**2))

        alpha1_roll = 2 * zeta_roll * wn_roll
        alpha2_roll = wn_roll**2
        alpha1_yaw = 2 * zeta_yaw * wn_yaw
        alpha2_yaw = wn_yaw**2

        self.kd_roll = alpha1_roll * P.J1x
        self.kp_roll = alpha2_roll * P.J1x
        self.kd_yaw = alpha1_yaw / b_yaw
        self.kp_yaw = alpha2_yaw / b_yaw

        print('kd_roll: ', self.kd_roll)
        print('kp_roll: ', self.kp_roll)

        print('kd_yaw: ', self.kd_yaw)
        print('kp_yaw: ', self.kp_yaw)
        print('ki_yaw: ', self.ki_yaw)
        # ***********************************************************

        # ******************* Discrete Time Calcs *******************
        # Sample Rate of Controller
        self.Ts = P.Ts

        # Dirty Derivative Values
        sigma = 0.01  # cutoff freq for dirty derivative
        self.beta = (2 * sigma - self.Ts) / (2 * sigma + self.Ts)

        # Delayed Variables
        # ----- Pitch (Theta) -----
        self.theta_d1 = 0
        self.theta_dot = 0
        self.error_theta_d1 = 0
        self.integrator_theta = 0

        # ----- Roll (Phi) -----
        self.phi_d1 = 0
        self.phi_dot = 0
        self.error_phi_d1 = 0

        # ----- Yaw (Psi) -----
        self.psi_d1 = 0
        self.psi_dot = 0
        self.error_psi_d1 = 0
        self.integrator_psi = 0

        # ***********************************************************

    def update(self, r, y):
        # **************** Pull Out Known States/Refs ****************
        phi = y[0][0]
        theta = y[1][0]
        psi = y[2][0]

        theta_ref = r[0][0]
        psi_ref = r[1][0]
        # ************************************************************

        # ********************* Pitch Controller *********************
        # Compute Errors
        error_theta = theta_ref - theta

        # Update Derivative/Integrator
        self.theta_dot = self.beta * self.theta_dot + (1 - self.beta) * ((theta - self.theta_d1) / P.Ts)
        self.integrator_theta = self.integrator_theta + (P.Ts / 2) * (error_theta + self.error_theta_d1)
        
        # Force Control
        force_unsat = (self.kp_pitch * error_theta) + (self.ki_pitch * self.integrator_theta) - (self.kd_pitch * self.theta_dot) + self.Fe
        force = saturate(force_unsat, -P.force_max, P.force_max)
        # ************************************************************

        # ********************* Roll/Yaw Controller *********************
        # ----- Outer Loop -----
        # Compute Errors
        error_psi = psi_ref - psi

        # Update Derivative/Integrator
        self.psi_dot = self.beta * self.psi_dot + (1 - self.beta) * ((psi - self.psi_d1) / P.Ts)
        self.integrator_psi = self.integrator_psi + (P.Ts / 2) * (error_psi + self.error_psi_d1)

        # Get Phi Ref
        phi_ref = (self.kp_yaw * error_psi) + (self.ki_yaw * self.integrator_psi) - (self.kd_yaw * self.psi_dot)

        # ----- Outer Loop -----
        # Compute Errors
        error_phi = phi_ref - phi

        # Update Derivative
        self.phi_dot = self.beta * self.phi_dot + (1 - self.beta) * ((phi - self.phi_d1) / P.Ts)

        # Get Torque
        torque_unsat = (self.kp_roll * error_phi) - (self.kd_roll * self.phi_dot)
        torque = saturate(torque_unsat, -P.torque_max, P.torque_max)
        # ***************************************************************

        # ********************* Return PWM and States & Update Variables *********************
        # Convert Force and Torque to PWM
        pwm = np.array([[force + torque / P.d],               # u_left
                      [force - torque / P.d]]) / (2 * P.km)   # r_right          
        pwm = saturate(pwm, 0, 1)

        # Update Delayed Variables
        self.theta_d1 = theta
        self.error_theta_d1 = error_theta
        self.psi_d1 = psi
        self.error_psi_d1 = error_psi
        self.phi_d1 = theta
        self.error_phi_d1 = error_phi

        # Return PWM and References
        return pwm, np.array([[phi_ref], [theta_ref], [psi_ref]])
        # ************************************************************************************


def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        if u > up_limit:
            u = up_limit
        if u < low_limit:
            u = low_limit
    else:
        for i in range(0, u.shape[0]):
            if u[i][0] > up_limit:
                u[i][0] = up_limit
            if u[i][0] < low_limit:
                u[i][0] = low_limit
    return u



