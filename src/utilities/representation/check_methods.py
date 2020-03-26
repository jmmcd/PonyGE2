from algorithm.parameters import params
from representation import individual
import numpy as np


def check_ind(ind, check):
    """
    Check all shallow aspects of an individual to ensure everything is correct.
    
    :param ind: An individual to be checked.
    :return: False if everything is ok, True if there is an issue.
    """

    if ind.genome == []:
        # Ensure all individuals at least have a genome.
        return True

    if ind.invalid and \
            ((check == "crossover" and params['NO_CROSSOVER_INVALIDS']) or
             (check == "mutation" and params['NO_MUTATION_INVALIDS'])):
        # We have an invalid.
        return True

    elif params['MAX_TREE_DEPTH'] and ind.depth > params['MAX_TREE_DEPTH']:
        # Tree is too deep.
        return True

    elif params['MAX_TREE_NODES'] and ind.nodes > params['MAX_TREE_NODES']:
        # Tree has too many nodes.
        return True

    elif params['MAX_GENOME_LENGTH'] and len(ind.genome) > \
            params['MAX_GENOME_LENGTH']:
        # Genome is too long.
        return True
