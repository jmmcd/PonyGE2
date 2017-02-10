from math import isnan

import numpy as np
np.seterr(all="raise")
import scipy

from algorithm.parameters import params
from utilities.fitness.get_data import get_data
from fitness.default_fitness import default_fitness
from utilities.fitness.math_functions import plog, rlog, pdiv, psqrt, aq, ppow, psqrt2, ppow2

class supervised_learning:
    """
    Fitness function for supervised learning, ie regression and
    classification problems. Given a set of training or test data,
    returns the error between y (true labels) and yhat (estimated
    labels).

    We can pass in the error metric and the dataset via the params
    dictionary. Of error metrics, MSE is suitable for regression,
    while F1-score, hinge-loss and others are suitable for
    classification.
    """

    maximise = False

    def __init__(self):
        # Get training and test data
        self.training_in, self.training_exp, self.test_in, self.test_exp = \
            get_data(params['DATASET'])

        # Find number of variables.
        self.n_vars = np.shape(self.test_in)[0]

        # Regression/classification-style problems use training and test data.
        self.training_test = True

    def __call__(self, ind, dist="training"):
        """
        Note that math functions used in the solutions are imported from either
        utilities.fitness.math_functions or called from numpy.

        :param ind: An individual to be evaluated.
        :param dist: An optional parameter for problems with training/test
        data. Specifies the distribution (i.e. training or test) upon which
        evaluation is to be performed.
        :return: The fitness of the evaluated individual.
        """

        phen = ind.phenotype

        if dist == "training":
            x = self.training_in
            y = self.training_exp
        elif dist == "test":
            x = self.test_in
            y = self.test_exp
        else:
            raise ValueError("Unknown dist: " + dist)

        try:
            if params['OPTIMIZE_CONSTANTS']:
                # if we are training, then optimize the constants by
                # gradient descent and save the resulting phenotype
                # string as ind.phenotype_with_c0123 (eg x[0] +
                # c[0] * x[1]**c[1]) and values for constants as
                # ind.opt_consts (eg (0.5, 0.7). Later, when testing,
                # use the saved string and constants to evaluate.
                if dist == "training":
                    fitness = self.optimize_constants(x, y, ind)
                else:
                    phen = ind.phenotype_with_c0123
                    c = ind.opt_consts
                    # phen will refer to x (ie test_in), and possibly to c
                    yhat = eval(phen)
                    assert np.isrealobj(yhat)

                    # let's always call the error function with the
                    # true values first, the estimate second
                    fitness = params['ERROR_METRIC'](y, yhat)

            else:
                # phenotype won't refer to C
                yhat = eval(phen)
                assert np.isrealobj(yhat)

                # let's always call the error function with the true
                # values first, the estimate second
                fitness = params['ERROR_METRIC'](y, yhat)


        except (FloatingPointError, ZeroDivisionError, OverflowError):
            # FP err can happen through eg overflow (lots of pow/exp calls)
            # ZeroDiv can happen when using unprotected operators
            fitness = default_fitness(self.maximise)
        except Exception as err:
            # other errors should not usually happen (unless we have
            # an unprotected operator) so user would prefer to see them
            print(err)
            raise

        # don't use "not fitness" here, because what if fitness = 0.0?!
        if isnan(fitness):
            fitness = default_fitness(self.maximise)

        return fitness


    def optimize_constants(self, x, y, ind):
        """Use gradient descent to search for values for the constants in
        ind.phenotype which minimise loss."""

        loss = params['ERROR_METRIC']
        s, n_consts = self.replace_ci_with_c0123(ind.phenotype)

        # let's save the old phenotype (has c[i] everywhere)
        ind.phenotype_original = ind.phenotype
        # and save the new one (has c[0], c[1], etc)
        ind.phenotype_with_c0123 = s

        f = eval("lambda x, c: " + s)

        if n_consts == 0:
            # ind doesn't refer to c: no need to optimize
            c = []
            fitness = loss(y, f(x, c))
            ind.opt_consts = c
            return fitness

        obj = lambda c: loss(y, f(x, c)) # obj is now a function of c
        # only for L-BFGS-B, using 0 as the init seems a reasonable
        # choice. But for scipy.curve_fit we might use [1.0] *
        # n_consts. Maybe other minimizers do better with some other
        # choices? There are other methods to try out.
        init = [0.0] * n_consts
        res = scipy.optimize.minimize(obj, init, method="L-BFGS-B")
        # the result is accessed like a dict
        ind.opt_consts = res['x'] # the optimum values of the constants

        # the most useful form of the phenotype: c[0], c[1] etc replaced
        # with actual numbers, so can be eval'd directly
        ind.phenotype = self.replace_c0123_with_values(s, ind.opt_consts)
        return res['fun'] # the value of the error metric at those values

    def replace_ci_with_c0123(self, s):
        # map occurences of c[i] to c[0], c[1], etc (and count number
        # of them). we always use c[i] in the grammar, because if we
        # used actual values like c[0], c[1], etc, the grammar could
        # create an individual with c[7] but no c[1].
        i = 0
        while "c[i]" in s:
            # replace the next occurrence only
            s = s.replace("c[i]", "c["+str(i)+"]", 1)
            i += 1
        n_consts = i
        # print("after c[i]: ", s)
        return s, n_consts

    def replace_c0123_with_values(self, s, c):
        for i in range(len(c)):
            s = s.replace("c[%d]" % i, str(c[i]))
        return s
