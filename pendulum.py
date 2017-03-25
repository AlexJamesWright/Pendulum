"""
Project
-------

To develop a code that can solve for the motion of a pendulum with the idea of
extending to a double pendulum and further. Aim is to investigate the effect of
initial conditions on the path.


---------------------------------------------
            \  ) ---> angle, theta
             \/
              \
               \
                \
                 \ }--> length, l
                  \
                   \
                    \
                     O
                    force, mg

theta.. = g cos(theta)
theta.(n+1) = theta.(n) + theta..(n) * dt
theta(n+1) = theta(n) + theta.(n) * dt


"""

import numpy as np
import pylab as plt

class simulation(object):
    def __init__(self, timeEv, initConds, length=1, mass=1, g=9.81):
        self.timeEv = timeEv            # Time evolution method
        self.initConds = initConds      # Initial data: (angle, angular velocity)
        self.length = length            # Length of pendulum (shaft is light)
        self.mass = mass                # Mass at end of pendulum
        self.coords = self.initConds    # Current (angle, angular velocity)
        self.allCoords = None             # [theta, thetadot, t] for every timestep
        self.g = g                      # Acceleration due to gravity
        self.t = 0                      # Time run
        self.dt = 0.001                 # Timestep - to be allocated dynamically
        self.iter = 0
        self.hi = 'Hiiiii'


    def runSim(self, endTime):
        """
        Run the simulation until endTime
        """
        self.allCoords = np.zeros((int(endTime/self.dt + 2), 3))
        self.allCoords[0] = self.coords[0], self.coords[1], self. t
        self.endTime = endTime
        while self.t < endTime:
            self.step(endTime)


    def step(self, endTime):
        if endTime - self.t < self.dt:
            self.dt = endTime - self.t
        self.coords = self.timeEv(self)
        self.t += self.dt
        self.allCoords[self.iter+1] = self.coords[0], self.coords[1], self.t
        self.iter += 1



    def plotAgainstTime(self):
        plt.figure()
        plt.plot(self.allCoords[:, 2], self.allCoords[:, 0])
        plt.show()


def euler(sim):
    theta, thetadot = sim.coords
    thetadot += sim.g * np.cos(theta)
    theta += thetadot * sim.dt
    return theta, thetadot



if __name__ == '__main__':
    initConds = np.array([-np.pi/2+1e-8, 0])
    sim = simulation(euler, initConds, 1, 1)
    sim.runSim(1.2)
    sim.plotAgainstTime()










