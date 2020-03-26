from PonyGE2.parameters import params
from PonyGE2.fitness import base_ff
from PonyGE2.data import arc_from_json


import itertools
import numpy as np

from PonyGE2.arc_utils import *


def matrix_dist(x, y):
    mx, nx = x.shape
    my, ny = y.shape
    pad_x = (max(0, my-mx), max(0, ny-nx))
    pad_y = (max(0, mx-my), max(0, nx-ny))
    
    norms = np.zeros((np.abs(mx-my)+1, np.abs(nx-ny)+1))
    for x_bottom, x_right in itertools.product(range(pad_x[0]+1), range(pad_x[1]+1)):
        for y_bottom, y_right in itertools.product(range(pad_y[0]+1), range(pad_y[1]+1)):
            x_ = np.pad(x, ((pad_x[0]-x_bottom, x_bottom), (pad_x[1]-x_right, x_right)))
            y_ = np.pad(y, ((pad_y[0]-y_bottom, y_bottom), (pad_y[1]-y_right, y_right)))
            print((x_bottom,x_right), (y_bottom,y_right))
            print(x_, y_)
            norms[max(x_bottom, y_bottom), max(x_right, y_right)] = np.linalg.norm(x_-y_, 'fro')
    
    return np.min(norms)
            

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
        self.x_train,
            self.y_train,
            self.x_test,
            self.y_test = arc_from_json(params['DATASET_TRAIN'])
        
        self.maximise = False

    def evaluate(self, ind, **kwargs):
        lambda_exp = 'lambda m, x, y: ' + ind.phenotype
        f = eval(lambda_exp)
        g = lambda x: apply(x, f)

        total_accuracy = sum(matrix_dist(g(x), y) for x, y in zip(self.x_train, self.y_train))
        return total_accuracy / float(len(self.x_train))
