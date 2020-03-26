from utilities import set_params
from ponyge import mane
from arc_evaluate import apply
import numpy as np

def pad(x):
    s = str(x)
    if x < 10:
        s = "0" + s
    if x < 100:
        s = "0" + s

    return s

def solve_task(task_nr):
    parameters = [
        "--cache", 
        "--codon_size", "1000", 
        "--crossover_probability", "0.5", 
        "--generations", "500", 
        "--max_genome_length", "500", 
        "--grammar_file", "arc.bnf", 
        "--initialisation", "operators.uniform_tree", 
        "--max_init_tree_depth", "5", 
        "--max_tree_depth", "30", 
        "--mutation_probability", "0.25", 
        "--population_size", "10",
        "--multicore",
        "--fitness_function", "arc_evaluate", 
        "--dataset_train", "arc/{}/Train.txt".format(pad(task_nr))
    ]

    set_params(parameters)

    individual, fitness = mane()
    individual = "lambda m, x, y: " + individual
    g = lambda x: apply(x, eval(individual))

    return g, individual, fitness

print(solve_task(52))