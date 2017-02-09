from algorithm.parameters import params

import editdistance # https://pypi.python.org/pypi/editdistance
import lzstring # https://pypi.python.org/pypi/lzstring/
import dtw # https://pypi.python.org/pypi/dtw

"""

This fitness function is for a sequence-match problem: we're given
an integer sequence target, say [0, 5, 0, 5, 0, 5], and we try to synthesize a program
(loops, if-statements, etc) which will *yield* that sequence, one item at a time.

There are several components of the fitness:

* concerning the program:
** length of the program (shorter is better)
** compressability of the program (non-compressible, ie DRY, is better)

* concerning distance from the target:
** dynamic time warping distance from the program's output to the target (lower is better)
** Levenshtein distance from the program's output to the target (lower is better).

"""

# available for use in synthesized programs
def succ(n, maxv=6): return min(n+1, maxv)
def pred(n, minv=0): return max(n-1, minv)

# the program will yield one item at a time, potentially forever.
# we only up to n items.
def truncate(n, g):
    for i in range(n):
        yield next(g)

# numerical difference, used as a component in DTW
def dist(t0, x0):
    return abs(t0 - x0)

# Dynamic time warping distance between two sequences
def dtw_dist(s, t):
    s = list(map(int, s))
    t = list(map(int, t))
    d, M, C, path = dtw.dtw(s, t, dist)
    return d

# Levenshtein distance between two sequences, normalised by length
# of the target -- hence this is *asymmetric*, not really a distance.
# don't normalise by length of the longer, because it would encourage
# evolution to create longer and longer sequences.
def lev_dist(s, t):
    return editdistance.eval(s, t) / len(s)

# convert to a string and compress. lzstring is a special-purpose compressor,
# more suitable for short strings than typical compressors.
def compress(s):
    s = ''.join(map(str, s))
    return lzstring.LZString().compress(s)

# compressability is in [0, 1]. It's high when the compressed string
# is much shorter than the original.
def compressability(s):
    return 1 - len(compress(s)) / len(s)

# program length is measured in characters, but in order to keep the values
# in a similar range to that of compressability, DTW and Levenshtein, we
# divide by 100. This is a bit arbitrary.
def proglen(s):
    return len(s) / 100.0

class sequence_match:
    def __init__(self):
        # --target will be a sequence such as (0, 5, 0, 5)
        self.target = eval(params['TARGET'])
        # we assume --extra_fitness_parameters is a comma-separated kv sequence, eg:
        # "alpha=0.5, beta=0.5, gamma=0.5"
        # which we can pass to the dict() constructor
        extra_fit_params = eval("dict("+params['EXTRA_FITNESS_PARAMETERS']+")")
        self.alpha = extra_fit_params['alpha']
        self.beta = extra_fit_params['beta']
        self.gamma = extra_fit_params['gamma']
        self.maximise = False

    def __call__(self, ind):
        # ind.phenotype will be a string incl fn defns etc. when we
        # exec it will create a value XXX_output_XXX, but we exec
        # inside an empty dict for safety.  but we put a couple of
        # useful primitives in the dict too.
        p = ind.phenotype
        d = {'pred': pred, 'succ': succ}
        exec(p, d)
        s = d['XXX_output_XXX'] # this is the program's output: a generator
        s = list(truncate(len(self.target), s)) # truncate the generator and convert to list
        t = self.target # our target


        # various weightings of four aspects of our fitness. the formula is:
        # fitness = gamma * dist + (1 - gamma) * length
        # where dist = alpha * lev_dist(t, s) + (1 - alpha) * dtw_dist(t, s)
        # and length = beta * proglen(t) + (1 - beta) * compressability(t)
        # but when any of alpha, beta and gamma is 0 or 1, we can save some calculation:

        if self.gamma > 0.0:
            if self.alpha > 0.0:
                lev_dist_v = lev_dist(t, s)
            else:
                lev_dist_v = 0.0
            if self.alpha < 1.0:
                dtw_dist_v = dtw_dist(t, s)
            else:
                dtw_dist_v = 0.0
            dist_v = self.alpha * lev_dist_v + (1 - self.alpha) * dtw_dist_v
        else:
            dist_v = 0.0

        if self.gamma < 1.0:
            if self.beta > 0.0:
                proglen_v = proglen(p)
            else:
                proglen_v = 0.0
            if self.beta < 1.0:
                compressability_v = compressability(p)
            else:
                compressability_v = 0.0
            length_v = self.beta * proglen_v + (1 - self.beta) * compressability_v
        else:
            length_v = 0.0

        return self.gamma * dist_v + (1 - self.gamma) * length_v


if __name__ == "__main__":
    # TODO write some tests here
    pass
