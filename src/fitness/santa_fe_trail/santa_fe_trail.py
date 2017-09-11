"""
Trail class, ant simulator and fitness function for santa fe trail
"""

from algorithm.parameters import params
from fitness.base_ff_classes.base_ff import base_ff
import copy
from os import path
from fitness.santa_fe_trail.gp import AntSimulator
#import pypeg2

"""
Ant simulator to simulate ants in Santa Fe trail environment
"""
    
class santa_fe_trail(base_ff):
    #Fitness function for Santa Fe Trail problem. Starts the AntSimulator
    #with the individuals phenotype and returns the amount of food eaten. 
    #Total food eaten is 89.

    def __init__(self):
        # Initialise base fitness function class.
        super().__init__()

        self.maximise = True

    def evaluate(self, ind, **kwargs):
        # Initialise the AntSimulator with maximum of 600 steps
        ant = AntSimulator(600)

        # Store the phenotype of the individual into code. The phenotype
        # python program sniplet that can be executed by AntSimulator
        code = ind.phenotype

        # Build python partial functions
        routine = ant.build_routine(code)

        # Pass the partial functions to the ant simulator
        ant.run(routine)

        # Return the amount of food eaten
        return ant.eaten