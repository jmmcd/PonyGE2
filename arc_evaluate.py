from PonyGE2.parameters import params
from PonyGE2.fitness import base_ff
from PonyGE2.data import arc_from_json

from numpy.random import normal
import numpy as np

from PonyGE2.arc_utils import *

def apply(x, f):
    result = np.zeros((x.shape[1], x.shape[0]))

    for row in range(x.shape[1]):
        for col in range(x.shape[0]):
            result[row, col] = f(x, row, col)

    return result

def accuracy(x, y, a=0.5):
    correct_shape_pixels = np.mean((x > 0) == (y > 0))
    correct_color_pixels = np.mean(x == y)

    return a * correct_shape_pixels + (1 - a) * correct_color_pixels

class arc_evaluate(base_ff):
    """Fitness function for matching a string. Takes a string and returns
    fitness. Penalises output that is not the same length as the target.
    Penalty given to individual string components which do not match ASCII
    value of target."""

    def __init__(self):
        # Initialise base fitness function class.
        super().__init__()
        self.x_train, self.y_train, self.x_test, self.y_test = arc_from_json(params['DATASET_TRAIN'])
       
        self.maximise = True

    def evaluate(self, ind, **kwargs):
        lambda_exp = 'lambda m, x, y: ' + ind.phenotype
        f = eval(lambda_exp)
        g = lambda x: apply(x, f)

        total_accuracy = sum(accuracy(g(x), y) for x, y in zip(self.x_train, self.y_train))
        return total_accuracy / float(len(self.x_train))