from sys import maxsize


def default_fitness(maximise):
    """ Return default fitness given maximization of minimization"""
    if maximise:
        return -maxsize
    else:
        return maxsize
