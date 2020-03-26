import numpy as np
from os import listdir
from os.path import join

from PonyGE2.utilities import init_params
from PonyGE2.ponyge import mane
from PonyGE2.arc_evaluate import apply
from PonyGE2.parameters import params

import os

def pad(x):
    s = str(x)
    if x < 10:
        s = "0" + s
    if x < 100:
        s = "0" + s

    return s

training_path = "training/"
evaluation_path = "evaluation/"
test_path = "test/"

training_tasks = sorted(listdir(training_path))
evaluation_tasks = sorted(listdir(evaluation_path))
test_tasks = sorted(listdir(test_path))

def __solve():
    init_params()

    individual, fitness = mane()
    individual = "lambda m, x, y: " + individual
    g = lambda x: apply(x, eval(individual))

    return g, individual, fitness

def solve_json(json_name):
    params["DATASET_TRAIN"] = join(training_path, json_name)
    params["DATASET_TEST"] = join(test_path, json_name)

    return __solve()

def solve_training_task(task_nr):

    params["DATASET_TRAIN"] = join(training_path, training_tasks[task_nr])
    params["DATASET_TEST"] = join(test_path, test_tasks[task_nr])

    return __solve()

#print(solve_training_task(52))