import operator
import random
import math
import unittest
import numpy as np

from deap import base
from deap import creator
from deap import tools

import math
import numpy as np
import pandas as pd
from copy import deepcopy
from deap import base, creator, algorithms, tools
import MGSurvE as srv
import unittest

def setup_stats(): 
    """
    Set up Statistics
    (note: we moved this outside because deap's tools did not work inside the class)
    """
    stats = tools.Statistics(lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)
    return stats

class Particle_Swarm:
    global trpMsk, immovable_trap_loc_mask
    trpMsk = None
    immovable_trap_loc_mask = None

    def __init__(
            self, traps, p_min, p_max, lnd, 
            num_particles=50, num_gens=500, 
            s_min=-3, s_max=3,
            phi1=2, phi2=2,
            optimFunction=srv.getDaysTillTrapped,
            optimFunctionArgs={'outer': np.mean, 'inner': np.max}
        ):
        """
        Initializes Particle Swarm Optimization. 
        
        Inputs: 
            traps: dataframe representing original locations of traps with the following 4 columns:
                x: x coordinate
                y: y coordinate
                t: type of trap
                f: 0 if movable, 1 if immovable 
            num_particles: number of particles
            num_gens: number of generations
            p_min: minimum longitude/latitude
            p_max: maximum longitude/latitude
            s_min: minimum speed of particle
            s_max: maximum speed of particle
            lnd: landscape defined by user
        """
        self.traps = traps
        self.num_particles = num_particles
        self.num_gens = num_gens  
        self.s_min = s_min
        self.s_max = s_max
        self.p_min = p_min
        self.p_max = p_max
        self.phi1 = phi1
        self.phi2 = phi2
        self.lnd = deepcopy(lnd)
        self.optimFunction = optimFunction
        self.optimFunctionArgs = optimFunctionArgs

        self.num_traps = traps.shape[0]
        Particle_Swarm.immovable_trap_loc_mask = self.getImmovableLocMask(traps)
        Particle_Swarm.trpMsk = srv.genFixedTrapsMask(self.lnd.trapsFixed, 2) * 1
        self.toolbox = self.setup_toolbox()

    def generate(size, p_min, p_max, smin, smax):
        """
        Generate 1 particle.

        Particle location is randomized for movable traps; it remains the same for immovable traps. 
        Particle speed is randomized for all elements in vector regardless of movable/immovable. This is fine because we multiply the speed by the trpMsk.
        smin, smax set by variables passed in.
        """
        random_locs = np.random.uniform(low=p_min, high=p_max, size=size) 
        random_locs_movable = np.multiply(random_locs, Particle_Swarm.trpMsk)
        random_with_original_immovable_locs = random_locs_movable + np.array(Particle_Swarm.immovable_trap_loc_mask)
        part = creator.Particle(list(random_with_original_immovable_locs))
        part.speed = [random.uniform(smin, smax) for _ in range(size)] # generate random SIZE speeds for all the particles above that are between SMIN and SMAX 
        part.smin = smin
        part.smax = smax
        return part

    def updateParticle(part, best, phi1, phi2):
        """
        Update location of the particle for 1 timestep.

        Inputs: 
            part: Array containing particle positions in (x, y) coordinates
            best: Global best particle position so far
            phi1: Upper bound for random number generation for particle velocity (u1)
            phi2: Upper bound for random number generation for particle velocity (u2)
            equation from https://machinelearningmastery.com/a-gentle-introduction-to-particle-swarm-optimization/
            speed/velocity per particle: V[t]
            X[t+1] = X[t] + V[t+1]
            V[t+1] = w * V[t] + c1 * r1 * (personal_best - X[t]) + c2 * r2 * (global_best - X[t])
                w: inertia weight constant = 1
                c1 = 1
                r1 = u1
                c2 = 1
                r2 = u2
                v_u1 = r1 * (personal_best - X[t]) = u1 * (part.best - part)
                v_u2 = r2 * (global_best - X[t]) = u2 * (best - part)
        """
        # calculate V[t + 1]
        u1 = (random.uniform(0, phi1) for _ in range(len(part)))
        u2 = (random.uniform(0, phi2) for _ in range(len(part)))
        v_u1 = map(operator.mul, u1, map(operator.sub, part.best, part))
        v_u2 = map(operator.mul, u2, map(operator.sub, best, part))
        part.speed = list(map(operator.add, part.speed, map(operator.add, v_u1, v_u2))) # update speed
        
        # if new speeds are outside speed limits, change to limit (keeping sign)
        for i, speed in enumerate(part.speed): 
            if abs(speed) < part.smin:
                part.speed[i] = math.copysign(part.smin, speed)
            elif abs(speed) > part.smax:
                part.speed[i] = math.copysign(part.smax, speed)

        ## WILL REMOVE WHEN FINISHED TESTING 
        # print(part)
        
        # update particles; calculate X[t+1]
        part[:] = list( map(operator.add, part, map(operator.mul, part.speed, Particle_Swarm.trpMsk)))
    
    
    def evaluate(self):
        """
        Main function that executes PSO.
        
        Output:
            pop: The resulting particle locations 
            logbook: Log of the changes that particles underwent
            best: The global best trap locations found by PSO
            
        """
        creator.create("FitnessMin", base.Fitness, weights=(-1.0, )) # want to minimize so that mosquitos die lol
        creator.create("Particle", list, fitness=creator.FitnessMin, speed=list, 
            smin=None, smax=None, best=None)
            
        pop = self.toolbox.population(n=self.num_particles)
        stats = setup_stats()

        logbook = tools.Logbook()
        logbook.header = ["gen", "evals"] + stats.fields
        best = None

        for g in range(self.num_gens):
            # update particle_best and global_best based off of evaluate function 
            for part in pop:
                part.fitness.values = self.toolbox.evaluate(part)
                if not part.best or part.best.fitness < part.fitness:
                    part.best = creator.Particle(part)
                    part.best.fitness.values = part.fitness.values
                if not best or best.fitness < part.fitness:
                    best = creator.Particle(part)
                    best.fitness.values = part.fitness.values
        
            # update particle's position 
            for part in pop:
                self.toolbox.update(part, best)

            # gather all the fitnesses in one list and print the stats
            logbook.record(gen=g, evals=len(pop), **stats.compile(pop))
            print(logbook.stream)
        log = pd.DataFrame(logbook)
        return pop, log, best

    def getImmovableLocMask(self, traps_df):
        """
        Generate a mask for immovable trap locations.
        If trap movable, location mask should be (0, 0).
        If trap immovable, location mask should be its original location (x, y)

        Input: 
            traps: dataframe representing original locations of traps with the following 4 columns:
                x: x coordinate
                y: y coordinate
                t: type of trap
                f: 0 if movable, 1 if immovable 
        """
        x = traps_df['x']
        y = traps_df['y']
        f = traps_df['f']
        x_y_locs = []
        for i in range(len(x)):
            if (f[i] == 0):
                x_y_locs.append([0, 0])
            else:
                x_y_locs.append([x[i], y[i]])
        x_y_locs = np.array(x_y_locs)
        return x_y_locs.flatten()

    def setup_toolbox(self):
        """
        Set up toolbox by registering particles, population, update function, and evaluate function. 

        Output: toolbox object for DEAP package 
        """
        toolbox = base.Toolbox()
        toolbox.register("particle", 
            Particle_Swarm.generate,size=2*self.num_traps, 
            p_min=self.p_min, p_max=self.p_max, 
            smin=self.s_min, smax=self.s_max
        )
        toolbox.register("population", 
            tools.initRepeat, list, toolbox.particle
        ) 
        toolbox.register("update", 
            Particle_Swarm.updateParticle, 
            phi1=self.phi1, phi2=self.phi2
        )
        toolbox.register("evaluate", 
            srv.calcFitness, 
            landscape=self.lnd,
            optimFunction=self.optimFunction,
            optimFunctionArgs=self.optimFunctionArgs
        )
        return toolbox


###############################################################################
# Main
###############################################################################
if __name__ == '__main__':
    unittest.main()