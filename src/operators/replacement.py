from algorithm.parameters import params


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


def steady_state(new_pop, old_pop):
    """
    Combines both old and new populations and returns the top
    POPULATION_SIZE individuals. In theory this is not true steady state
    replacement, but it is functionally equivalent.
     
    :param new_pop: The new population (e.g. after selection, variation, &
    evaluation).
    :param old_pop: The previous generation population.
    :return: The 'POPULATION_SIZE' new population.
    """
    
    # Combine both populations
    total_pop = old_pop + new_pop
    
    # Sort the combined population
    total_pop.sort(reverse=True)
    
    # Return the top POPULATION_SIZE individuals of the combined population.
    return total_pop[:params['POPULATION_SIZE']]
