from collections import deque
from os import path

import numpy as np
from multiprocessing import Pool

from PonyGE2.representation import Tree
from PonyGE2.parameters import params
from PonyGE2.fitness import evaluate_fitness
from PonyGE2.operators import crossover, mutation, replacement, steady_state, selection, initialisation
from PonyGE2.stats import stats, get_stats
from PonyGE2.python_filter import python_filter
from PonyGE2.utilities import pool_init

"""
step.py
"""

def step(individuals):
    """
    Runs a single generation of the evolutionary algorithm process:
        Selection
        Variation
        Evaluation
        Replacement
    
    :param individuals: The current generation, upon which a single
    evolutionary generation will be imposed.
    :return: The next generation of the population.
    """

    # Select parents from the original population.
    parents = selection(individuals)

    # Crossover parents and add to the new population.
    cross_pop = crossover(parents)

    # Mutate the new population.
    new_pop = mutation(cross_pop)

    # Evaluate the fitness of the new population.
    new_pop = evaluate_fitness(new_pop)

    # Replace the old population with the new population.
    individuals = replacement(new_pop, individuals)

    # Generate statistics for run so far
    get_stats(individuals)
    
    return individuals


def steady_state_step(individuals):
    """
    Runs a single generation of the evolutionary algorithm process,
    using steady state replacement.

    :param individuals: The current generation, upon which a single
    evolutionary generation will be imposed.
    :return: The next generation of the population.
    """
    
    individuals = steady_state(individuals)
    
    return individuals

"""
search_loop.py
"""

def search_loop():
    """
    This is a standard search process for an evolutionary algorithm. Loop over
    a given number of generations.
    
    :return: The final population after the evolutionary process has run for
    the specified number of generations.
    """

    if params['MULTICORE']:
        # initialize pool once, if mutlicore is enabled
        params['POOL'] = Pool(processes=params['CORES'], initializer=pool_init,
                              initargs=(params,))  # , maxtasksperchild=1)

    # Initialise population
    individuals = initialisation(params['POPULATION_SIZE'])

    # Evaluate initial population
    individuals = evaluate_fitness(individuals)

    # Generate statistics for run so far
    get_stats(individuals)

    # Traditional GE
    for generation in range(1, (params['GENERATIONS']+1)):
        stats['gen'] = generation

        # New generation
        individuals = params['STEP'](individuals)

    if params['MULTICORE']:
        # Close the workers pool (otherwise they'll live on forever).
        params['POOL'].close()

    return individuals