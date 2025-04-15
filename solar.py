import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
from numpy.linalg import norm
import seaborn as sns; sns.set()
from planet import Planet
from sklearn.metrics import r2_score
class Solar(object):
    """
    Class to run the orbital simulation
    """

    def __init__(self):
        """
        The initialization of the Sun, the planets, and the system
        """
        inputdata = []
        filein = open("/Users/janghyoin/Desktop/2024-1/Computer Simulation/parameters-solar.txt", "r", encoding = "UTF-8")
        for line in filein.readlines():
            if not line.startswith("#"):
                inputdata.append(line)
        filein.close()

        # simulation parameters
        self.niter = int(inputdata[0])
        self.dt = float(inputdata[1])
        self.G = float(inputdata[2])

        # list for mars and moons
        self.bodies = []
        
        # rest of input data is mars and moon data in four line "chunks"
        # first entry must be mars
        for i in range(3, len(inputdata)-4, 4):
            name = inputdata[i]
            mass = float(inputdata[i+1])
            orbit = float(inputdata[i+2])
            colour = inputdata[i+3]
            self.bodies.append(Planet(name, mass, orbit, colour))
                                
        # set initial positions and velocities relative to sun
        # sun must be first element in bodies list!
        for i in range(0, len(self.bodies)):
            self.bodies[i].initialise(self.G, self.bodies[0])
        
        # dictionary to store the name of the planets and their orbital periods
        self.orbital_periods = {}

        # lists for planetary alignments
        self.alignment_detected = []
        self.mean_angle = []

        # create an array for patches (planet and moons)
        self.patches = []
    
    def init(self):
        # initialiser for animator
        return self.patches

    def animate(self, i):
        """
        Update the position and velocity at the _th time step. 
        Check the new year, the total energy, and the planetary alignment.
        """
        # keep track of time in earth years
        time = (i+1)*self.dt

        # update positions
        for j in range(0, len(self.bodies)):
            self.bodies[j].updatePos(self.G, self.dt)
            self.patches[j].center = self.bodies[j].r
            
        # then update velocities
        for j in range(0, len(self.bodies)):
            for k in range(0, len(self.bodies)):
                if j != k:
                    self.bodies[j].updateVel(self.G, self.dt, self.bodies[k])

        # check year and print year if new year for any planet (except the Sun)
        for j in range(1, len(self.bodies)):
            if self.bodies[j].newYear(self.bodies[0]):
                print (f"{self.bodies[j].name.strip()} "
                       f"{self.bodies[j].year} years = {time} earth years.")

                # if new year is earth year, also print total energy
                if self.bodies[j].name.strip() == "earth":
                    # need to convert from earth masses AU^2 yr^-2 to kg m^2 s-2 (J)
                    # 1 earth mass = 5.97219e24 kg
                    # 1 AU = 1.496e+11 m
                    c =(5.97219e+24*1.496e+11*1.496e+11)/(3.154e+7*3.154e+7)
                    energy = self.energy()*c
                    print(f"Time = {time} earth years. "
                          f"Total energy = {energy:.3e} J.")    

        # detect the planetary alignment and print time and the mean angle of the planets if the planets align      
        alignment, mean_angle = self.planetary_alignment()
        if alignment:
            print(f"Planetary alignment occured at time {time:.3f} earth years with mean angle {mean_angle:.3f}")
        
        return self.patches

    def runSimulation(self, i):
        """
        Simulate but does not print out anything.
        Check the new year and planetary alignments.
        """
        # keep track of time in earth years
        time = (i+1)*self.dt

        # update positions
        for j in range(0, len(self.bodies)):
            self.bodies[j].updatePos(self.G, self.dt)
            
        # then update velocities
        for j in range(0, len(self.bodies)):
            for k in range(0, len(self.bodies)):
                if j != k:
                    self.bodies[j].updateVel(self.G, self.dt, self.bodies[k])

        # check year and print year if new year for any planet (except the Sun)
        for j in range(1, len(self.bodies)):
            if self.bodies[j].newYear(self.bodies[0]):
                planet_name = self.bodies[j].name.strip()
                # Save the orbital period information in the dictionary {planet name:orbital periods} at the first new year
                if planet_name in self.orbital_periods:
                    continue
                self.orbital_periods[planet_name] = time                   
        # check the planetary alignments
        alignment, mean_angle = self.planetary_alignment()
        if alignment:
            self.alignment_detected.append(time) # check the time of the simulation
            self.mean_angle.append(mean_angle) # if the planets align, append the mean angle to the list mean_angle. Otherwise, append 0.
        else:
            self.alignment_detected.append(time)
            self.mean_angle.append(0)
    
    def energy(self):
        """
        Calculate the total energy in the summation of the kinetic and the potential energy.
        """
        ke = 0.0 
        pe = 0.0
        for j in range(0, len(self.bodies)):
            ke += self.bodies[j].kineticEnergy()
            for k in range(0, len(self.bodies)):
                if k != j:
                    r = norm(self.bodies[k].r - self.bodies[j].r)
                    pe -= self.G*self.bodies[j].m*self.bodies[k].m / r
        # divide pe by two to avoid double countin
        pe = pe / 2
        totEnergy = ke + pe
        return totEnergy

    def calcTotalEnergy(self, i):
        """
        Calculate the total energy in the summation of the kinetic and the potential energy.
        """
        ke = 0.0
        pe = 0.0
        for j in range(0, len(self.bodies)):
            ke += self.bodies[j].kineticEnergy()
            for k in range(0, len(self.bodies)):
                if k != j:
                    r = norm(self.bodies[k].r - self.bodies[j].r)
                    pe -= self.G*self.bodies[j].m*self.bodies[k].m / r
        # divide pe by two to avoid double countin
        pe = pe / 2
        totEnergy = ke + pe
        print(f"Time = {i} iterations. Total energy = {totEnergy:.3e} J") 

    def run(self):
        # set up the plot components        
        fig = plt.figure()
        ax = plt.axes()


        # get orbital radius of outermost moon to set size of
        # orbiting bodies and of plot
        # hacky - should really check to see which moon is outermost

        maxOrb = math.sqrt(np.dot(self.bodies[-1].r, self.bodies[-1].r))

        # add the planet and moons to the Axes and patches
        for i in range(0, len(self.bodies)):
            if i == 0:
                self.patches.append(
                    ax.add_patch(plt.Circle(self.bodies[i].r, 0.05*maxOrb,
                                            color=self.bodies[i].c, animated=True)))
            else:
                self.patches.append(
                    ax.add_patch(plt.Circle(self.bodies[i].r, 0.02*maxOrb,
                                            color=self.bodies[i].c, animated=True)))

        # set up the axes
        # scale axes so circle looks like a circle and set limits
        # with border b for prettier plot

        b = 1.2
        lim = maxOrb*b
        print(lim)
        ax.axis("scaled")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)

        anim = FuncAnimation(
            fig, self.animate, init_func=self.init,
            frames=self.niter, repeat=False, interval=1, blit=True)

        plt.show()

    def beeman_energy(self, years):
        """
        Calculate total energy of the orbital system using Beeman method.
        """
        # Convert years to the number of iterations
        iterations = int(years / self.dt)
        c =(5.97219e+24*1.496e+11*1.496e+11)/(3.154e+7*3.154e+7)        
        # Run the simulation with the Beeman method
        energies = []
        for i in range(iterations):
            self.runSimulation(i) # update position and velocities using runSimulation method
            # Calculate total energy with Beeman method
            total_energy = self.energy() * c
            energies.append(total_energy)
        # Return calculated energy list
        return energies
    
    def plotEnergyConservation(self, years):
        """
        Plot the total energy vs years of the orbital system using Beeman method.
        """
        iterations = int(years / self.dt)
        x = np.linspace(0, years, num = iterations)
        plt.plot(x, self.beeman_energy(years))
        plt.xlabel('Time (Earth years)')
        plt.ylabel('Total Energy')
        plt.title('Energy vs. Time using Beeman Methods')
        plt.show()         

    def planetary_alignment(self, threshold_degrees=5):
        """
        Detect planetary alignment based on mean angle of angles between planets.
        """
        # Compute pairwise angles between vectors pointing from the Sun to different planets
        angles = []
        for i in range(len(self.bodies)):
            for j in range(i + 1, len(self.bodies)): # range from (i + 1) to avoid calculating the angle twice
                vector1 = self.bodies[i].r - self.bodies[0].r  # Vector pointing from the Sun to planet i
                vector2 = self.bodies[j].r - self.bodies[0].r  # Vector pointing from the Sun to planet j
                if norm(vector1) == 0 or norm(vector2) == 0: # to avoid RuntimeWarning: invalid value encountered in scalar divide error
                    angle_radians = 0  
                else:
                    # use arccosine and convert it to radian to obtain the angles between the planets
                    angle_radians = np.arccos(np.dot(vector1, vector2) / (norm(vector1) * norm(vector2)))
                angles.append(np.degrees(angle_radians))

        # Calculate the mean angle
        mean_angle = np.mean(angles)

        # Check if the planets align within +-5 degrees of the mean angle
        # If the planets are aligned, the method returns 'True' and the mean angle of the angles
        alignment_detected = all(abs(a - mean_angle) < threshold_degrees for a in angles)
        return alignment_detected, mean_angle

    def plotAlignmentYears(self, years):
        """
        Plot the years when planetary alignments occur.
        """
        iterations = int(years / self.dt) # Convert years to the number of iterations
        for i in range(iterations): # run simulation
            self.runSimulation(i)
        
        # print the year and the mean angle of that time if the planets align
        for i in range(len(self.alignment_detected)):
            if self.mean_angle[i] != 0 :
                print(f"The planetary alignment occurs at {self.alignment_detected[i]} Earth years with the mean angle {self.mean_angle[i]:.3f}")
        
        # draw the plot using year as the x-axis, the mean angle as the y-axis
        # when the planets align, the mean angle is indicated on the plot
        plt.plot(self.alignment_detected,self.mean_angle, '.', label='Alignment Year')
        plt.xlabel('Earth Years')
        plt.ylabel('Mean angles')
        plt.title('Planetary Alignment Years')
        plt.grid(True)
        plt.show()

    def plotOrbitalPeriods(self, years):
        """
        Plot the orbital periods of each planet.
        """
        # Convert years to the number of iterations
        iterations = int(years / self.dt)   

        # Run the simulation with the Beeman method
        for i in range(iterations):
            self.runSimulation(i)    

        # Create a list of actual orbits of the planets in year (Mercury, Venus, Earth, Mars, Jupiter)
        # Reference: Royal Museums Greenwich
        self.actual_orbital_periods = [87.97/365.20, 224.70/365.20, 365.20/365.20, 686.98/365.20, 11.86]
        
        # Get planet names and corresponding orbital periods
        planet_names = list(self.orbital_periods.keys())
        orbital_periods = [self.orbital_periods[i] for i in planet_names]

        width = 0.2
        n_planet = np.arange(len(planet_names))

        # draw two bar graphs (one for the simulated orbital period, one for the actual orbital period) side by sdie
        plt.figure(1)                     
        plt.bar(n_planet, orbital_periods, width, label='Simulated orbital periods')
        plt.bar(n_planet + width, self.actual_orbital_periods, width, label = 'Actual orbital periods')
        plt.xlabel('Planets')
        plt.ylabel('Orbital Period (Earth years)')
        plt.title('Orbital Period of Each Planet')
        plt.xticks(n_planet + width / 2, planet_names, rotation=30)
        plt.legend(loc='best')
        
        # draw the scattered plot setting simulated orbital period as the x-axis and the actual period as the y-axis
        # impose on the regression line
        plt.figure(2)
        periods = np.array([min(orbital_periods), max(orbital_periods)])
        # get the regression line
        fit_line = np.polyfit(orbital_periods, self.actual_orbital_periods, 1)
        fit_y = periods * fit_line[0] + fit_line[1]
        # get the R^2
        est_y = np.array(orbital_periods) * fit_line[0] + fit_line[1] 
        r2 = r2_score(self.actual_orbital_periods, est_y)
        plt.scatter(orbital_periods, self.actual_orbital_periods, color = 'r')
        plt.plot(periods, fit_y, color = 'orange')
        plt.text(6, 5, '$R^2$ = %.4f'%r2, size = 10)
        plt.text(6, 4, 'y = %.4fx + %d'%(fit_line[0], fit_line[1]), size = 10)
        plt.xlabel('Simulated orbital period')
        plt.ylabel('Actual orbital period')
        plt.title('Orbital Periods')
        
        plt.show()

class SolarEuler(Solar):
    """
    Class to run the Euler method. Inherited from the 'Solar' class.
    """
    def __init__(self):
        Solar.__init__(self)     

    def euler_energy(self, years):
        # Convert years to iterations based on the simulation parameters
        iterations = int(years / self.dt)
        c =(5.97219e+24*1.496e+11*1.496e+11)/(3.154e+7*3.154e+7)
        
        # Run the simulation with the Euler method
        energies = []
        for i in range(iterations):
            # Update positions and velocities using Euler method
            for body in self.bodies:
                body.updatePos_euler(self.dt)

            for body in self.bodies:
                for other_body in self.bodies:
                    if body != other_body:
                        body.updateVel_euler(self.G, self.dt, other_body)
            # Calculate total energy with Euler method
            total_energy = self.energy() * c
            energies.append(total_energy)
        
        return energies
        
    def euler_energy_symplecticeuler(self, years): #the symplectic Euler method
        # Convert years to iterations 
        iterations = int(years / self.dt)
        c =(5.97219e+24*1.496e+11*1.496e+11)/(3.154e+7*3.154e+7)
        energies = []
        for i in range(iterations):
            # Update positions and velocities using the symplectic Euler method
            # Update position first (same to the Direct Euler)
            for body in self.bodies:
                body.updatePos_euler(self.dt)
            # Update velocity
            for body in self.bodies:
                for other_body in self.bodies:
                    if body != other_body:
                        body.updateVel_symplecticeuler(self.G, self.dt, other_body)
            # Calculate total energy with Euler method
            total_energy = self.energy() * c
            energies.append(total_energy)
        
        return energies
    
    def energy_comparison(self, years):
        iterations = int(years / self.dt)
        solar = Solar() # state two different objects of each methods to run the simulation at the same time
        solar_Euler = SolarEuler()
        beeman_energies = solar.beeman_energy(years) # calculate the total energy using the Beeman methods
        euler_energies = solar_Euler.euler_energy(years) # calculate the total energy using the Direct Euler methods

        x = np.linspace(0, years, num = iterations)
        plt.plot(x, euler_energies, label='Euler')
        plt.plot(x, beeman_energies, label='Beeman')
        plt.xlabel('Time (Earth years)')
        plt.ylabel('Total Energy')
        plt.title('Comparison of Total Energy over Time (Beeman vs. Euler)')
        plt.legend()
        plt.show() 

    def energy_comparison_symplecticeuler(self, years):
        iterations = int(years / self.dt)
        solar = Solar() # state two different objects of each methods to run the simulation at the same time
        solar_Euler = SolarEuler()
        beeman_energies = solar.beeman_energy(years) # calculate the total energy using the Beeman methods
        euler_energies = solar_Euler.euler_energy_symplecticeuler(years) # calculate the total energy using the symplectic Euler methods

        x = np.linspace(0, years, num = iterations)
        plt.plot(x, euler_energies, label='symplectic Euler')
        plt.plot(x, beeman_energies, label='Beeman')
        plt.xlabel('Time (Earth years)')
        plt.ylabel('Total Energy')
        plt.title('Comparison of Total Energy over Time (Beeman vs. symplectic Euler)')
        plt.legend()
        plt.show() 

if __name__ == "__main__":
    # Create Solar object
    solar = Solar()
    # Run the simulation
    solar.run()
    # Exp 1
    # Check the orbital periods of the planets. Run the simulation for 20 Earth years.
    solar.plotOrbitalPeriods(20)
    # Exp 2
    # Check the total energy of the orbital system (Beeman). Run the simulation for 800 Earth years.
    solar.plotEnergyConservation(800)
    # Create SolarEuler object
    euler = SolarEuler()
    # Compare the total energy of the orbital systems (Beeman vs. Direct Euler). Run the simulation for 800 Earth years.    
    euler.energy_comparison(800)
    # Compare the total energy of the orbital systems (Beeman vs. symplectic Euler). Run the simulation for 800 Earth years.
    euler.energy_comparison_symplecticeuler(800)
    # Exp 4
    # Check the planetary alignments. Run the simulation for 10,000 years.
    solar.plotAlignmentYears(10000)