from algorithm.mapper import mapper, genome_map, map_tree_from_genome
from fitness.fitness import default_fitness
from algorithm.parameters import params
from random import randint


class Individual(object):
    """A GE individual"""

    def __init__(self, genome, ind_tree):

        self.phenotype, self.genome, self.tree, self.nodes, self.invalid, \
        self.depth, self.used_codons = mapper(genome, ind_tree)
        self.fitness = default_fitness(params['FITNESS_FUNCTION'].maximise)
        self.name = None

    def __lt__(self, other):
        if params['FITNESS_FUNCTION'].maximise:
            return self.fitness < other.fitness
        else:
            return other.fitness < self.fitness

    def __str__(self):
        return ("Individual: " +
                str(self.phenotype) + "; " + str(self.fitness))

    def evaluate(self, dist="training"):
        """ Evaluates phenotype in fitness function on either training or test
        distributions and sets fitness"""

        if params['PROBLEM'] in ("regression", "classification"):
            # The problem is regression, e.g. has training and test data
            self.fitness = params['FITNESS_FUNCTION'](self.phenotype, dist)
        else:
            self.fitness = params['FITNESS_FUNCTION'](self.phenotype)

        if params['MULTICORE']:
            return self
