import numpy as np 
import massParam as P


class massDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],      # initial angle
            [P.zdot0]    # initial angular rate
        ])  
        # Mass of the arm, kg
        self.m = P.m * (1.+alpha*(2.*np.random.rand()-1.))
        # Damping coefficient, Ns
        self.b = P.b * (1.+alpha*(2.*np.random.rand()-1.))  
        # the spring constant, N/m
        self.k = P.k * (1.+alpha*(2.*np.random.rand()-1.))
        # sample rate at which the dynamics are propagated
        self.Ts = P.Ts  
        self.force_limit = P.F_max

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
        zd = state[1][0]
        zdd = (F - self.b*zd - self.k*z)/self.m
        xdot = np.array([[zd], [zdd]])
        return xdot

    def h(self):
        # return the output equations
        # could also use input u if needed
        z = self.state[0][0]
        y = np.array([[z]])
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
