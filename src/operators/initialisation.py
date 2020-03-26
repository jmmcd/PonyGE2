from math import floor
from os import path, getcwd, listdir
from random import shuffle, randint

from algorithm.parameters import params
from representation import individual
from representation.derivation import generate_tree
from representation.individual import Individual
from representation.tree import Tree
from utilities.representation.python_filter import python_filter


def initialisation(size):
    """
    Perform selection on a population in order to select a population of
    individuals for variation.
    
    :param size: The size of the required population.
    :return: A full population generated using the specified initialisation
    technique.
    """

    # Decrease initialised population size by the number of seed individuals
    # (if any) to ensure that the total initial population size does not exceed
    # the limit.
    size -= len(params['SEED_INDIVIDUALS'])

    # Initialise empty population.
    individuals = params['INITIALISATION'](size)

    # Add seed individuals (if any) to current population.
    individuals.extend(params['SEED_INDIVIDUALS'])

    return individuals
    

def uniform_tree(size):
    """
    Create a population of individuals by generating random derivation trees.
     
    :param size: The size of the required population.
    :return: A full population composed of randomly generated individuals.
    """
    
    return [generate_ind_tree(params['MAX_TREE_DEPTH'],
                              "random") for _ in range(size)]


def generate_ind_tree(max_depth, method):
    """
    Generate an individual using a given subtree initialisation method.

    :param max_depth: The maximum depth for the initialised subtree.
    :param method: The method of subtree initialisation required.
    :return: A fully built individual.
    """

    # Initialise an instance of the tree class
    ind_tree = Tree(str(params['BNF_GRAMMAR'].start_rule["symbol"]), None)

    # Generate a tree
    genome, output, nodes, _, depth = generate_tree(ind_tree, [], [], method,
                                                    0, 0, 0, max_depth)

    # Get remaining individual information
    phenotype, invalid, used_cod = "".join(output), False, len(genome)

    if params['BNF_GRAMMAR'].python_mode:
        # Grammar contains python code

        phenotype = python_filter(phenotype)

    # Initialise individual
    ind = individual.Individual(genome, ind_tree, map_ind=False)

    # Set individual parameters
    ind.phenotype, ind.nodes = phenotype, nodes
    ind.depth, ind.used_codons, ind.invalid = depth, used_cod, invalid

    # Generate random tail for genome.
    ind.genome = genome + [randint(0, params['CODON_SIZE']) for
                           _ in range(int(ind.used_codons / 2))]

    return ind

