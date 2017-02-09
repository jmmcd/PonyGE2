from algorithm.parameters import params


class string_match:
    """Fitness function for matching a string. Takes a string and returns
    fitness. Penalises output that is not the same length as the target.
    Penalty given to individual string components which do not match ASCII
    value of target."""

    maximise = False

    def __init__(self):
        # Set target string.
        self.target = params['TARGET']

    def __call__(self, ind):
        guess = ind.phenotype
        fitness = max(len(self.target), len(guess))
        # Loops as long as the shorter of two strings
        for (t_p, g_p) in zip(self.target, guess):
            if t_p == g_p:
                # Perfect match.
                fitness -= 1
            else:
                # Imperfect match, find ASCII distance to match.
                fitness -= 1 / (1 + (abs(ord(t_p) - ord(g_p))))
        return fitness
