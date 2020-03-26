from fitness.evaluation import evaluate_fitness
from algorithm.parameters import params
from operators.crossover import crossover_inds
from operators.selection import selection

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

