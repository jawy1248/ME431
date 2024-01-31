import numpy as np 
import blockbeamParam as P


class blockbeamDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],         # initial position
            [P.theta0],     # initial angle
            [P.zdot0],      # initial velocity
            [P.thetadot0]   # initial angular velocity
        ])  
        # Mass of the arm, kg
        self.m1 = P.m1 * (1.+alpha*(2.*np.random.rand()-1.))
        self.m2 = P.m2 * (1.+alpha*(2.*np.random.rand()-1.))
        self.g = P.g
        self.ell = P.length
        # sample rate at which the dynamics are propagated
        self.Ts = P.Ts  
        self.force_limit = P.Fmax

    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input torque
        u = saturate(u, self.force_limit)
        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y

    def f(self, state, u):
        # Return xdot = f(x,u), the system state update equations
        # re-label states for readability
        F = u

        z = state[0][0]
        theta = state[1][0]
        zd = state[2][0]
        thetad = state[3][0]

        zdd = -self.g*np.sin(theta) + z*(thetad**2)
        thetadd = (6*F*self.ell*np.cos(theta) - 3*self.ell*self.g*self.m2*np.cos(theta) - 6*self.g*self.m1*z*np.cos(theta) - 12*self.m1*z*zd*thetad) / (2*(self.ell**2)*self.m2 + 6*self.m1*(z**2))

        xdot = np.array([[zd], [thetad], [zdd], [thetadd]])
        return xdot

    def h(self):
        # return the output equations
        # could also use input u if needed
        z = self.state[0][0]
        theta = self.state[1][0]
        y = np.array([[z], [theta]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = (self.f(self.state, u)).astype(float)
        F2 = (self.f(self.state + self.Ts / 2 * F1, u)).astype(float)
        F3 = (self.f(self.state + self.Ts / 2 * F2, u)).astype(float)
        F4 = (self.f(self.state + self.Ts * F3, u)).astype(float)
        self.state = self.state + (self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4))

    
def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
