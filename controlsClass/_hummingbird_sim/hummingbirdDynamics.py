import numpy as np 
import hummingbirdParam as P


class HummingbirdDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.phi0],  # roll angle
            [P.theta0],  # pitch angle
            [P.psi0],  # yaw angle
            [P.phidot0],  # roll rate
            [P.thetadot0],  # pitch rate
            [P.psidot0],  # yaw rate
        ])

        # vary the actual physical parameters
        self.ell1 = P.ell1 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ell2 = P.ell2 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ell3x = P.ell3x * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ell3y = P.ell3y * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ell3z = P.ell3z * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ellT = P.ellT * (1. + alpha * (2. * np.random.rand() - 1.))
        self.d = P.d * (1. + alpha * (2. * np.random.rand() - 1.))
        self.m1 = P.m1 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.m2 = P.m2 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.m3 = P.m3 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J1x = P.J1x * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J1y = P.J1y * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J1z = P.J1z * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J2x = P.J2x * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J2y = P.J2y * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J2z = P.J2z * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J3x = P.J3x * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J3y = P.J3y * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J3z = P.J3z * (1. + alpha * (2. * np.random.rand() - 1.))
 
    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input force
        u = saturate(u, P.torque_max)
        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y

    def f(self, state, u):
        # Return xdot = f(x,u)
        phi = state[0][0]
        theta = state[1][0]
        psi = state[2][0]
        phidot = state[3][0]
        thetadot = state[4][0]
        psidot = state[5][0]
        
        pwm_left = u[0][0]
        pwm_right = u[1][0]

        # For easy acces
        sinP = np.sin(phi)
        sinP2 = np.sin(phi)**2
        sinT = np.sin(theta)
        sinT2 = np.sin(theta)**2
        cosP = np.cos(phi)
        cosP2 = np.cos(phi)**2
        cosT = np.cos(theta)
        cosT2 = np.cos(theta)**2

        ml1 = self.m1*(self.ell1**2)
        ml2 = self.m2*(self.ell2**2)

        # The equations of motion go here
        M22 = ml1 + ml2 + self.J2y + (self.J1y*cosP2) + (self.J1z*sinP2)
        M23 = (self.J1y - self.J1z)*sinP*cosP*cosT
        M33 = ((ml1 + ml2 + self.J2z + (self.J1y*sinP2) + (self.J1z*cosP2))*cosT2) + ((self.J1x + self.J2x)*sinT2) + (self.m3*((self.ell3x**2) + (self.ell3y**2))) + self.J3z
        M = np.array([[self.J1x,           0,  -(self.J1x*sinT)],
                      [0,                 M22,       M23],
                      [-(self.J1x*sinT),  M23,       M33]
                      ])
        
        N33 = 2*(self.J1x + self.J2x - ml1 - ml2 - self.J2z - (self.J1y*sinP2) - (self.J1z*cosP2))*sinT*cosT

        line1 = (thetadot**2)*(self.J1z - self.J1y)*sinP*cosP*sinT
        line2 = (((self.J1y - self.J1z)*(cosP2 - sinP2)) - self.J1x)*cosT*phidot*thetadot
        line3 = ((self.J1z - self.J1y)*sinP*cosP*sinT*(thetadot**2)) + (2*(self.J1y - self.J1z)*sinP*cosP*phidot*psidot)
        line4 = 2*(-ml1 - ml2 - self.J2z + self.J1x + self.J2x + (self.J1y*sinP2) + (self.J1z*sinP2))*sinT*cosT*thetadot*psidot

        C = np.array([[(self.J1y - self.J1z)*sinP*cosP*((thetadot**2) - cosT2*(psidot**2)) + ((self.J1y - self.J1z)*(cosP2 - sinP2) - self.J1x)*cosT*thetadot*psidot],
                      [(2*(self.J1z - self.J1y)*sinP*cosP*phidot*thetadot) + ((((self.J1y - self.J1z)*(cosP2 - sinP2)) + self.J1x)*cosT*phidot*psidot) - (0.5*N33*(psidot**2))],
                      [line1 + line2 + line3 + line4],
                     ])
        
        partialP = np.array([[0],
                             [((self.m1*self.ell1) + (self.m2*self.ell2))*P.g*cosT],
                             [0],
                            ])
        
        force = P.km * (pwm_left + pwm_right)
        torque = self.d * P.km * (pwm_left - pwm_right)
        tau = np.array([[torque],
                        [self.ellT*(force)*cosP],
                        [(self.ellT*(force)*cosT*sinP) - (torque*sinT)]])
        
        beta = 0.001
        B = beta*np.eye(3)

        qddot = np.linalg.inv(M) @ (-C - partialP + tau - B @ state[3:6])
        phiddot = qddot[0][0]
        thetaddot = qddot[1][0]
        psiddot = qddot[2][0]
        # build xdot and return
        xdot = np.array([[phidot],
                         [thetadot],
                         [psidot],
                         [phiddot],
                         [thetaddot],
                         [psiddot]])
        return xdot

    def h(self):
        # return y = h(x)
        phi = self.state[0][0]
        theta = self.state[1][0]
        psi = self.state[2][0]
        y = np.array([[phi], [theta], [psi]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = (self.f(self.state, u)).astype(float)
        F2 = (self.f(self.state + P.Ts / 2 * F1, u)).astype(float)
        F3 = (self.f(self.state + P.Ts / 2 * F2, u)).astype(float)
        F4 = (self.f(self.state + P.Ts * F3, u)).astype(float)
        self.state = self.state + (P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4))

def saturate(u, limit):
    for i in range(0, u.shape[0]):
        if abs(u[i][0]) > limit:
            u[i][0] = limit * np.sign(u[i][0])
    return u
