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
from matplotlib import animation
from IPython.core.display import HTML

class simulation(object):
    def __init__(self, timeEv, initConds, length=np.array([[1]]), mass=np.array([[1]]), g=9.81, Npend=1):
        assert(initConds.shape[0] == Npend), 'Npends={}, but only {} initConds given'.format(Npend, initConds.shape[0])
        assert(length.shape[0] == Npend), 'Npends={}, but only {} lengths given'.format(Npend, length.shape[0])
        assert(mass.shape[0] == Npend), 'Npends={}, but only {} masses given'.format(Npend, mass.shape[0])
        self.g = g                      # Acceleration due to gravity
        self.Npends = Npend
        self.timeEv = timeEv            # Time evolution method
        self.initConds = initConds      # Initial data: (angle, angular velocity)
        self.length = length            # Length of pendulums (array-like)
        self.mass = mass                # Mass at end of pendulums (array-like)
        self.coords = np.zeros((Npend, 3))
        self.coords[:, :2] = self.initConds    # Current (angle, angular velocity)
        self.setInit() #set initial accelerations
        self.allCoords = np.zeros((Npend, 1, 4))         # [theta, thetadot, thetadotdot, t] for every timestep
        for i in range(Npend):
            self.allCoords[i, 0, :] = initConds[i, 0], initConds[i, 1], self.coords[i, 2], 0
        self.t = 0                      # Time run
        self.dtInit = 0.01
        self.dt = self.dtInit                 # Timestep - to be allocated dynamically
        self.iter = 0

    def setInit(self):
        for i in range(self.Npends):
            thetaSum = np.sum(self.coords[:i, 0])
            thetadotSum = np.sum(self.coords[:i, 1])
            self.coords[i, 2] = self.g * self.length[i] * np.cos(thetaSum) - thetadotSum

    def runSim(self, endTime):
        """
        Run the simulation until endTime
        """

        self.endTime = endTime
        while self.t < endTime:
            self.step(endTime)


    def step(self, endTime):
        if endTime - self.t < self.dt:
            self.dt = endTime - self.t
        thetaVec, thetadotVec, thetaddotVec = self.timeEv(self)
        for i in range(self.Npends):
            self.coords[i, :] = np.array([thetaVec[i], thetadotVec[i], thetaddotVec[i]])
        self.t += self.dt
        self.store(self.coords)
        self.iter += 1

    def store(self, coords):
        temp = np.zeros((self.Npends, self.allCoords.shape[1]+1, 4))
        temp[:, :-1, :] = self.allCoords
        temp[:, -1, :3] = self.coords[:, :]
        temp[:, -1, -1] = self.t
        self.allCoords = temp

    def plotAgainstTime(self, pend=-1):
        plt.figure()
        plt.xlim([-0.05, self.t])
        plt.plot(self.allCoords[pend, :, -1], self.allCoords[pend, :, 0])
        plt.show()

    def plotPendulums(self):
        plt.figure()
        maxwidth = self.length.sum() * 1.05
        plt.xlim([-maxwidth, maxwidth])
        plt.ylim([-maxwidth, maxwidth])
        plt.plot([-maxwidth, maxwidth], [0, 0], linewidth=5)
        posX = np.zeros(self.Npends + 1)
        posY = np.zeros(self.Npends + 1)
        for i in range(self.Npends):
            thetaSum = np.sum(self.allCoords[:i+1, -1, 0])
            posX[i+1] = posX[i] + self.length[i] * np.cos(thetaSum)
            posY[i+1] = posY[i] - self.length[i] * np.sin(thetaSum)
            if i < self.Npends:
                plt.plot([posX[i], posX[i+1]], [posY[i], posY[i+1]])
        plt.show()

    def animatePends(self):
        fig, ax = plt.subplots()
        maxwidth = np.sum(sim.length)*1.05
        ax.set_xlim(-maxwidth, +maxwidth)
        ax.set_ylim(-maxwidth, +maxwidth)
        line, = ax.plot([], [], '-0')
        time_template = 'time = %.1fs'
        time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

        def init():
            line.set_data([], [])
            time_text.set_text('')
            return line, time_text

        timeSteps = sim.allCoords[0, :, 0].shape[0]
        posX = np.zeros((timeSteps, sim.Npends + 1))
        posY = np.zeros((timeSteps,sim.Npends + 1))
        for k in range(sim.Npends):
            for j in range(timeSteps):
                thetaSum = np.sum(sim.allCoords[:k+1, j, 0])
                posX[j, k+1] = posX[j, k] + sim.length[k] * np.cos(thetaSum)
                posY[j, k+1] = posY[j, k] - sim.length[k] * np.sin(thetaSum)



        def animate(i):

            line.set_data([posX[i]], [posY[i]])
            time_text.set_text(time_template%(i*self.dtInit))
            return line, time_text

        anim = animation.FuncAnimation(fig, animate, init_func=init, frames=timeSteps, interval=self.dtInit*1000, blit=True)
        return HTML(anim.to_html5_video())



def euler(sim):
    thetaVec = np.zeros(Npends)
    thetadotVec = np.zeros(Npends)
    thetaddotVec = np.zeros(Npends)
    for i in range(sim.Npends):
        theta, thetadot, thetaddot = sim.coords[i]
        thetaSum = np.sum(sim.coords[:i+1, 0])
        thetaddotSum = np.sum(sim.coords[:i, 2])
        thetaddot = sim.length[i] * sim.g * np.cos(thetaSum) - thetaddotSum
        # And the reverse torque from the previous step
        if i < sim.Npends-1:
            COMfrac = sim.mass[i+1] * sim.length[i+1] / (sim.mass[i] + sim.mass[i+1])
            T1 = sim.mass[i+1] * sim.g * sim.length[i+1] * sim.coords[i+1, 2]
            T2 = T1 * COMfrac / (sim.length[i+1] - COMfrac)
            revTorq = np.cos(sim.coords[i+1, 0]) * T2
            thetaddot -= revTorq*sim.dt
        # And the second reverse torque on the pendulum below
        if i > 0:
            COMfrac = sim.mass[i] * sim.length[i] / (sim.mass[i] + sim.mass[i-1])
            Ttheta=sim.mass[i-1] * sim.length[i-1] * sim.g * np.cos(sim.coords[i-1, 0])
            T1 = np.cos(sim.coords[i, 0]) * Ttheta
            T2 = T1 * (sim.length[i] - COMfrac) / COMfrac
            revTorq = T2
            thetaddot -= revTorq

        thetadot += thetaddot*sim.dt
        theta += thetadot * sim.dt
        thetaVec[i] = theta
        thetadotVec[i] = thetadot
        thetaddotVec[i] = thetaddot
    return thetaVec, thetadotVec, thetaddotVec



#if __name__ == '__main__':

#Double pendulum
Npends=2
initConds = np.array([[np.pi/2, 0], [0.4, 0]])
length = np.array([[1], [1]])
mass = np.array([[1], [1]])
sim = simulation(euler, initConds, length, mass, Npend=Npends)
sim.runSim(20)
sim.animatePends()
#    sim.plotPendulums()


#for i in range(Npends):
#    sim.plotAgainstTime(i)











