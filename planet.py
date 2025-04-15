import math
import matplotlib
import numpy as np
from numpy.linalg import norm

class Planet(object):
    """
    Planet class
    """

    def __init__(self, name, mass, orbit, colour):
        self.name = name
        # mass in kg
        self.m = mass
        # orbital radius in m
        self.orbit = orbit
        # colour - need to strip trailing line return!
        self.c = colour.strip()
        # set year to zero
        self.year = 0

    def initialise(self, G, p):
        # inital position, initial coords = (orbit radius, 0)
        self.r = np.array([self.orbit, 0])
        # inital velocity, tangential to position
        # speed = sqrt(G*marsmass/r)
        if self.orbit == 0.0:
            self.v = np.array([0.0, 0])
        else:
            vel = math.sqrt(G*p.m/self.orbit)
            self.v = np.array([0.0, vel])
        # intial accelatation, using gravitational force law
        if self.orbit == 0.0:
            self.a = np.array([0, 0])
        else:
            self.a = self.updateAcc(G, p)
        # set acc_old = acc to start Beeman
        self.a_old = self.a

    def updatePos(self, G, dt):
        # keep old position to check for year
        self.r_old = self.r
        
        # update position first: Beeman
        self.r = self.r + self.v*dt + (4*self.a - self.a_old)*dt*dt/6.0
        
    def updateVel(self, G, dt, p):
        # update velocity second: Beeman
        a_new = self.updateAcc(G, p)
        self.v = self.v + + (2*a_new + 5*self.a - self.a_old)*dt/6.0
        # now update acc ready for next iteration
        self.a_old = self.a
        self.a = a_new

    def updateAcc(self, G, p):
        # update acc (gravitational force law)
        pos = self.r - p.r
        a = -G*p.m*pos/math.pow(norm(pos),3)
        return a
    
    def updatePos_euler(self, dt):
        # update the position first (the Direct Euler, the symplectic Euler)
        self.a_old = self.a # keep old acceleration for the velocity
        self.r += self.v * dt  
    
    def updateVel_euler(self, G, dt, p):
        # then, update the velocity and the acceleration (the Direct Euler)
        self.v += self.a_old * dt 
        self.a = self.updateAcc(G, p)
    
    def updateVel_symplecticeuler(self, G, dt, p):
        # update the velocity with new accelartaion (the symplectic Euler)
        self.v += self.updateAcc(G, p) * dt

    def newYear(self, Sun):
        # update the year when the planet passes the + x axis of Sun
        if self.r_old[1] < Sun.r[1] and self.r[1] >= Sun.r[1]:
            self.year +=1
            return True
        else:
            return False

    # determine kinetic energy
    def kineticEnergy(self):
        # ke in J
        ke = (np.dot(self.v, self.v))*self.m/2
        return ke