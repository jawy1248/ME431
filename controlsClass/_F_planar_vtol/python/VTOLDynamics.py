import numpy as np 
import VTOLParam as P

class Dynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.z0],         # initial distance
            [P.h0],         # initial height
            [P.theta0],     # initial angle
            [P.zdot0],      # initial velocity right
            [P.hdot0],      # initial velocity up
            [P.thetadot0]   # initial angular velocity
        ])  
        # Params from VTOL Params
        self.mc = P.mc * (1.+alpha*(2.*np.random.rand()-1.))
        self.mr = P.mr * (1.+alpha*(2.*np.random.rand()-1.))
        self.Jc = P.Jc
        self.d = P.d
        self.mu = P.mu
        self.g = P.g
        # sample rate at which the dynamics are propagated
        self.Ts = P.Ts  
        self.force_limit = P.fmax

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
        fr = u[0][0]
        fl = u[1][0]

        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        zd = state[3][0]
        hd = state[4][0]
        thetad = state[5][0]

        zdd = (-np.sin(theta)*(fr+fl) - self.mu*zd + P.F_wind)/(2*self.mr + self.mc)
        hdd = (np.cos(theta)*(fr+fl) - self.g*(2*self.mr + self.mc))/(2*self.mr + self.mc)
        thetadd = (self.d*(fr-fl))/(self.Jc + 2*(self.d**2)*self.mr)

        xdot = np.array([[zd], [hd], [thetad], [zdd], [hdd], [thetadd]])
        return xdot

    def h(self):
        # return the output equations
        # could also use input u if needed
        z = self.state[0][0]
        h = self.state[1][0]
        theta = self.state[2][0]
        y = np.array([[z], [h], [theta]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = (self.f(self.state, u)).astype(float)
        F2 = (self.f(self.state + self.Ts / 2 * F1, u)).astype(float)
        F3 = (self.f(self.state + self.Ts / 2 * F2, u)).astype(float)
        F4 = (self.f(self.state + self.Ts * F3, u)).astype(float)
        self.state = self.state + (self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4))

    
def saturate(u, limit):
    if abs(u[0][0]) > limit:
        u[0][0] = limit * np.sign(u[0][0])
    if abs(u[1][0]) > limit:
        u[1][0] = limit * np.sign(u[1][0])
    return u
