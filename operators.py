from math import floor
from os import getcwd, listdir, path
from random import choice, randint, random, sample, shuffle

from PonyGE2.parameters import params
from PonyGE2.fitness import evaluate_fitness
from PonyGE2.representation import generate_tree
from PonyGE2.representation import Individual
from PonyGE2.representation import Tree
from PonyGE2.python_filter import python_filter


"""
initilisation.py
"""

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
    ind = Individual(genome, ind_tree, map_ind=False)

    # Set individual parameters
    ind.phenotype, ind.nodes = phenotype, nodes
    ind.depth, ind.used_codons, ind.invalid = depth, used_cod, invalid

    # Generate random tail for genome.
    ind.genome = genome + [randint(0, params['CODON_SIZE']) for
                           _ in range(int(ind.used_codons / 2))]

    return ind

"""
crossover.py
"""

def crossover(parents):
    """
    Perform crossover on a population of individuals. The size of the crossover
    population is defined as params['GENERATION_SIZE'] rather than params[
    'POPULATION_SIZE']. This saves on wasted evaluations and prevents search
    from evaluating too many individuals.
    
    :param parents: A population of parent individuals on which crossover is to
    be performed.
    :return: A population of fully crossed over individuals.
    """

    # Initialise an empty population.
    cross_pop = []
    
    while len(cross_pop) < params['GENERATION_SIZE']:
        
        # Randomly choose two parents from the parent population.
        inds_in = sample(parents, 2)

        # Perform crossover on chosen parents.
        inds_out = crossover_inds(inds_in[0], inds_in[1])
        
        if inds_out is None:
            # Crossover failed.
            pass
        
        else:
                        
            # Extend the new population.
            cross_pop.extend(inds_out)

    return cross_pop

def crossover_inds(parent_0, parent_1):
    """
    Perform crossover on two selected individuals.
    
    :param parent_0: Parent 0 selected for crossover.
    :param parent_1: Parent 1 selected for crossover.
    :return: Two crossed-over individuals.
    """

    # Create copies of the original parents. This is necessary as the
    # original parents remain in the parent population and changes will
    # affect the originals unless they are cloned.
    ind_0 = parent_0.deep_copy()
    ind_1 = parent_1.deep_copy()

    # Crossover cannot be performed on invalid individuals.
    if not params['INVALID_SELECTION'] and (ind_0.invalid or ind_1.invalid):
        s = "operators.crossover.crossover\nError: invalid individuals " \
            "selected for crossover."
        raise Exception(s)

    # Perform crossover on ind_0 and ind_1.
    inds = params['CROSSOVER'](ind_0, ind_1)

    # Check each individual is ok (i.e. does not violate specified limits).
    checks = [check_ind(ind, "crossover") for ind in inds]

    if any(checks):
        # An individual violates a limit.
        return None

    else:
        # Crossover was successful, return crossed-over individuals.
        return inds

def variable_onepoint(p_0, p_1):
    """
    Given two individuals, create two children using one-point crossover and
    return them. A different point is selected on each genome for crossover
    to occur. Note that this allows for genomes to grow or shrink in
    size. Crossover points are selected within the used portion of the
    genome by default (i.e. crossover does not occur in the tail of the
    individual).
    
    :param p_0: Parent 0
    :param p_1: Parent 1
    :return: A list of crossed-over individuals.
    """

    # Get the chromosomes.
    genome_0, genome_1 = p_0.genome, p_1.genome

    # Uniformly generate crossover points.
    max_p_0, max_p_1 = get_max_genome_index(p_0, p_1)
        
    # Select unique points on each genome for crossover to occur.
    pt_0, pt_1 = randint(1, max_p_0), randint(1, max_p_1)

    # Make new chromosomes by crossover: these slices perform copies.
    if random() < params['CROSSOVER_PROBABILITY']:
        c_0 = genome_0[:pt_0] + genome_1[pt_1:]
        c_1 = genome_1[:pt_1] + genome_0[pt_0:]
    else:
        c_0, c_1 = genome_0[:], genome_1[:]

    # Put the new chromosomes into new individuals.
    ind_0 = Individual(c_0, None)
    ind_1 = Individual(c_1, None)

    return [ind_0, ind_1]

def get_max_genome_index(ind_0, ind_1):
    """
    Given two individuals, return the maximum index on each genome across
    which operations are to be performed. This can be either the used
    portion of the genome or the entire length of the genome.
    
    :param ind_0: Individual 0.
    :param ind_1: Individual 1.
    :return: The maximum index on each genome across which operations are to be
             performed.
    """

    # Get length of entire genome.
    return  len(ind_0.genome), len(ind_1.genome)

# Set attributes for all operators to define linear or subtree representations.
variable_onepoint.representation = "linear"

"""
mutation.py
"""

def mutation(pop):
    """
    Perform mutation on a population of individuals. Calls mutation operator as
    specified in params dictionary.

    :param pop: A population of individuals to be mutated.
    :return: A fully mutated population.
    """

    # Initialise empty pop for mutated individuals.
    new_pop = []

    # Iterate over entire population.
    for ind in pop:

        # If individual has no genome, default to subtree mutation.
        if not ind.genome and params['NO_MUTATION_INVALIDS']:
            new_ind = subtree(ind)

        else:
            # Perform mutation.
            new_ind = params['MUTATION'](ind)

        # Check ind does not violate specified limits.
        check = check_ind(new_ind, "mutation")

        while check:
            # Perform mutation until the individual passes all tests.

            # If individual has no genome, default to subtree mutation.
            if not ind.genome and params['NO_MUTATION_INVALIDS']:
                new_ind = subtree(ind)

            else:
                # Perform mutation.
                new_ind = params['MUTATION'](ind)

            # Check ind does not violate specified limits.
            check = check_ind(new_ind, "mutation")

        # Append mutated individual to population.
        new_pop.append(new_ind)

    return new_pop


def int_flip_per_codon(ind):
    """
    Mutate the genome of an individual by randomly choosing a new int with
    probability p_mut. Works per-codon. Mutation is performed over the
    effective length (i.e. within used codons, not tails) by default;
    within_used=False switches this off.

    :param ind: An individual to be mutated.
    :return: A mutated individual.
    """

    # Set effective genome length over which mutation will be performed.
    eff_length = get_effective_length(ind)

    if not eff_length:
        # Linear mutation cannot be performed on this individual.
        return ind

    # Set mutation probability. Default is 1 over the length of the genome.
    if params['MUTATION_PROBABILITY'] and params['MUTATION_EVENTS'] == 1:
        p_mut = params['MUTATION_PROBABILITY']
    elif params['MUTATION_PROBABILITY'] and params['MUTATION_EVENTS'] > 1:
        s = "operators.mutation.int_flip_per_codon\n" \
            "Error: mutually exclusive parameters for 'MUTATION_PROBABILITY'" \
            "and 'MUTATION_EVENTS' have been explicitly set.\n" \
            "       Only one of these parameters can be used at a time with" \
            "int_flip_per_codon mutation."
        raise Exception(s)
    else:
        # Default mutation events per individual is 1. Raising this number
        # will influence the mutation probability for each codon.
        p_mut = params['MUTATION_EVENTS']/eff_length

    # Mutation probability works per-codon over the portion of the
    # genome as defined by the within_used flag.
    for i in range(eff_length):
        if random() < p_mut:
            ind.genome[i] = randint(0, params['CODON_SIZE'])

    # Re-build a new individual with the newly mutated genetic information.
    new_ind = Individual(ind.genome, None)

    return new_ind


def subtree(ind):
    """
    Mutate the individual by replacing a randomly selected subtree with a
    new randomly generated subtree. Guaranteed one event per individual, unless
    params['MUTATION_EVENTS'] is specified as a higher number.

    :param ind: An individual to be mutated.
    :return: A mutated individual.
    """

    def subtree_mutate(ind_tree):
        """
        Creates a list of all nodes and picks one node at random to mutate.
        Because we have a list of all nodes, we can (but currently don't)
        choose what kind of nodes to mutate on. Handy.

        :param ind_tree: The full tree of an individual.
        :return: The full mutated tree and the associated genome.
        """

        # Find the list of nodes we can mutate from.
        targets = ind_tree.get_target_nodes([], target=params[
                                          'BNF_GRAMMAR'].non_terminals)

        # Pick a node.
        new_tree = choice(targets)

        # Set the depth limits for the new subtree.
        if params['MAX_TREE_DEPTH']:
            # Set the limit to the tree depth.
            max_depth = params['MAX_TREE_DEPTH'] - new_tree.depth

        else:
            # There is no limit to tree depth.
            max_depth = None

        # Mutate a new subtree.
        generate_tree(new_tree, [], [], "random", 0, 0, 0, max_depth)

        return ind_tree

    if ind.invalid:
        # The individual is invalid.
        tail = []

    else:
        # Save the tail of the genome.
        tail = ind.genome[ind.used_codons:]

    # Allows for multiple mutation events should that be desired.
    for i in range(params['MUTATION_EVENTS']):
        ind.tree = subtree_mutate(ind.tree)

    # Re-build a new individual with the newly mutated genetic information.
    ind = Individual(None, ind.tree)

    # Add in the previous tail.
    ind.genome = ind.genome + tail

    return ind


def get_effective_length(ind):
    """
    Return the effective length of the genome for linear mutation.

    :param ind: An individual.
    :return: The effective length of the genome.
    """

    if not ind.genome:
        # The individual does not have a genome; linear mutation cannot be
        # performed.
        return None

    elif ind.invalid:
        # Individual is invalid.
        eff_length = len(ind.genome)

    elif params['WITHIN_USED']:
        eff_length = min(len(ind.genome), ind.used_codons)

    else:
        eff_length = len(ind.genome)

    return eff_length


int_flip_per_codon.representation = "linear"
subtree.representation = "subtree"

"""
replacement.py
"""

def replacement(new_pop, old_pop):
    """
    Given a new population and an old population, performs replacement using
    specified replacement operator.
    
    :param new_pop: Newly generated population (after selection, variation &
    evaluation).
    :param old_pop: The previous generation population.
    :return: Replaced population.
    """
    return params['REPLACEMENT'](new_pop, old_pop)

def generational(new_pop, old_pop):
    """
    Replaces the old population with the new population. The ELITE_SIZE best
    individuals from the previous population are appended to new pop regardless
    of whether or not they are better than the worst individuals in new pop.
    
    :param new_pop: The new population (e.g. after selection, variation, &
    evaluation).
    :param old_pop: The previous generation population, from which elites
    are taken.
    :return: The 'POPULATION_SIZE' new population with elites.
    """

    # Sort both populations.
    old_pop.sort(reverse=True)
    new_pop.sort(reverse=True)
    
    # Append the best ELITE_SIZE individuals from the old population to the
    # new population.
    for ind in old_pop[:params['ELITE_SIZE']]:
        new_pop.insert(0, ind)
    
    # Return the top POPULATION_SIZE individuals of the new pop, including
    # elites.
    return new_pop[:params['POPULATION_SIZE']]

def steady_state(individuals):
    """
    Runs a single generation of the evolutionary algorithm process,
    using steady state replacement:
        Selection
        Variation
        Evaluation
        Replacement
        
    Steady state replacement uses the Genitor model (Whitley, 1989) whereby
    new individuals directly replace the worst individuals in the population
    regardless of whether or not the new individuals are fitter than those
    they replace. Note that traditional GP crossover generates only 1 child,
    whereas linear GE crossover (and thus all crossover functions used in
    PonyGE) generates 2 children from 2 parents. Thus, we use a deletion
    strategy of 2.

    :param individuals: The current generation, upon which a single
    evolutionary generation will be imposed.
    :return: The next generation of the population.
    """

    # Initialise counter for new individuals.
    ind_counter = 0

    while ind_counter < params['POPULATION_SIZE']:
        
        # Select parents from the original population.
        parents = selection(individuals)

        # Perform crossover on selected parents.
        cross_pop = crossover_inds(parents[0], parents[1])
        
        if cross_pop is None:
            # Crossover failed.
            pass

        else:
            # Mutate the new population.
            new_pop = mutation(cross_pop)
        
            # Evaluate the fitness of the new population.
            new_pop = evaluate_fitness(new_pop)
    
            # Sort the original population
            individuals.sort(reverse=True)
    
            # Combine both populations
            total_pop = individuals[:-len(new_pop)] + new_pop
        
            # Increment the ind counter
            ind_counter += params['GENERATION_SIZE']

    # Return the combined population.
    return total_pop

"""
selection.py
"""

def selection(population):
    """
    Perform selection on a population in order to select a population of
    individuals for variation.

    :param population: input population
    :return: selected population
    """

    return params['SELECTION'](population)

def tournament(population):
    """
    Given an entire population, draw <tournament_size> competitors randomly and
    return the best. Only valid individuals can be selected for tournaments.

    :param population: A population from which to select individuals.
    :return: A population of the winners from tournaments.
    """

    # Initialise list of tournament winners.
    winners = []

    # The flag "INVALID_SELECTION" allows for selection of invalid individuals.
    if params['INVALID_SELECTION']:
        available = population
    else:
        available = [i for i in population if not i.invalid]

    while len(winners) < params['GENERATION_SIZE']:
        # Randomly choose TOURNAMENT_SIZE competitors from the given
        # population. Allows for re-sampling of individuals.
        competitors = sample(available, params['TOURNAMENT_SIZE'])

        # Return the single best competitor.
        winners.append(max(competitors))

    # Return the population of tournament winners.
    return winners

"""
check_methods.py
"""

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

def load_population(target):
    """
    Given a target folder, read all files in the folder and load/parse
    solutions found in each file.
    
    :param target: A target folder stored in the "seeds" folder.
    :return: A list of all parsed individuals stored in the target folder.
    """

    # Set path for seeds folder
    path_1 = path.join(getcwd(), "..", "seeds")

    if not path.isdir(path_1):
        # Seeds folder does not exist.
    
        s = "scripts.seed_PonyGE2.load_population\n" \
            "Error: `seeds` folder does not exist in root directory."
        raise Exception(s)
    
    path_2 = path.join(path_1, target)

    if not path.isdir(path_2):
        # Target folder does not exist.
    
        s = "scripts.seed_PonyGE2.load_population\n" \
            "Error: target folder " + target + \
            " does not exist in seeds directory."
        raise Exception(s)
    
    # Get list of all target individuals in the target folder.
    target_inds = [i for i in listdir(path_2) if i.endswith(".txt")]
       
    # Initialize empty list for seed individuals.
    seed_inds = []

    for ind in target_inds:
        # Loop over all target individuals.

        # Get full file path.
        file_name = path.join(path_2, ind)

        # Initialise None data for ind info.
        genotype, phenotype = None, None

        # Open file.
        with open(file_name, "r") as f:
            
            # Read file.
            raw_content = f.read()
            
            # Read file.
            content = raw_content.split("\n")
            
            # Check if genotype is already saved in file.
            if "Genotype:" in content:
                
                # Get index location of genotype.
                gen_idx = content.index("Genotype:") + 1
                
                # Get the genotype.
                try:
                    genotype = eval(content[gen_idx])
                except:
                    s = "scripts.seed_PonyGE2.load_population\n" \
                        "Error: Genotype from file " + file_name + \
                        " not recognized: " + content[gen_idx]
                    raise Exception(s)
            
            # Check if phenotype (target string) is already saved in file.
            if "Phenotype:" in content:
    
                # Get index location of genotype.
                phen_idx = content.index("Phenotype:") + 1
    
                # Get the phenotype.
                phenotype = content[phen_idx]
                
                # TODO: Current phenotype is read in as single-line only. Split is performed on "\n", meaning phenotypes that span multiple lines will not be parsed correctly. This must be fixed in later editions.
            
            elif "Genotype:" not in content:
                # There is no explicit genotype or phenotype in the target
                # file, read in entire file as phenotype.
                phenotype = raw_content

        if genotype:
            # Generate individual from genome.
            ind = Individual(genotype, None)
            
            if phenotype and ind.phenotype != phenotype:
                s = "scripts.seed_PonyGE2.load_population\n" \
                    "Error: Specified genotype from file " + file_name + \
                    " doesn't map to same phenotype. Check the specified " \
                    "grammar to ensure all is correct: " + \
                    params['GRAMMAR_FILE']
                raise Exception(s)
            
        # Add new ind to the list of seed individuals.
        seed_inds.append(ind)

    return seed_inds