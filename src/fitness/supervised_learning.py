from math import isnan

import numpy as np
np.seterr(all="raise")

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

        if dist == "test":
            x = self.test_in
            y = self.test_exp
        elif dist == "training":
            x = self.training_in
            y = self.training_exp

        try:
            yhat = eval(phen)  # phen will refer to "x", created above
            assert np.isrealobj(yhat)

            # let's always call the error function with the true values first,
            # the estimate second
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
