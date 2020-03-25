from algorithm.parameters import params
from fitness.base_ff_classes.base_ff import base_ff
from utilities.fitness.get_data import get_data

from numpy.random import normal
import numpy as np

from utilities.arc_utils import *

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
        self.x_train, self.y_train, self.x_test, self.y_test = get_data(params['DATASET_TRAIN'], params['DATASET_TEST'])
        self.x_train = self.x_train.T
       
        # Set target string.
        # self.target = params['TARGET']

    def evaluate(self, ind, **kwargs):
        lambda_exp = 'lambda m, x, y: ' + ind.phenotype
        f = eval(lambda_exp)
        g = lambda x: apply(x, f)

        inputs = []
        outputs = []
        for dp in self.x_train:
            x_len = int(dp[0] * dp[1])
            y_len = int(dp[2] * dp[3])
            reshaped_x = np.array(dp[4:(4 + x_len)]).reshape(int(dp[0]), int(dp[1]))
            reshaped_y = np.array(dp[(4 + y_len):]).reshape(int(dp[2]), int(dp[3]))
            inputs.append(reshaped_x)
            outputs.append(reshaped_y)

        total_accuracy = 0
        for x, y in zip(inputs, outputs):
            pred = g(x)
            total_accuracy += accuracy(pred, y) 

        print(total_accuracy / float(len(inputs)))
        return total_accuracy / float(len(inputs))
