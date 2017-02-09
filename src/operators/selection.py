from random import sample

from algorithm.parameters import params


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
        competitors.sort(reverse=True)
        winners.append(competitors[0])
    
    # Return the population of tournament winners.
    return winners


def truncation(population):
    """
    Given an entire population, return the best <proportion> of them.
    
    :param population: A population from which to select individuals.
    :return: The best <proportion> of the given population.
    """
    
    # Sort the original population.
    population.sort(reverse=True)
    
    # Find the cutoff point for truncation.
    cutoff = int(len(population) * float(params['SELECTION_PROPORTION']))
    
    # Return the best <proportion> of the given population.
    return population[:cutoff]
