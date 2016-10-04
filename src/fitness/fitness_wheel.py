from fitness.string_match import StringMatch
from fitness.regression import Regression
from fitness.progsys import ProgSys
from utilities.error_metrics import inverse_f1_score


def set_fitness_params(problem, params):
    if problem in ("regression", "classification"):
        return "grammars/" + params['SUITE'] + ".bnf", params['SUITE']
    elif problem == "string_match":
        return "grammars/letter.bnf", params['STRING_MATCH_TARGET']
    elif problem == "progsys":
        return "grammars/" + params['SUITE'] + ".bnf", params['SUITE']
    else:
        print("Error: Problem not specified correctly")
        exit(2)


def set_fitness_function(problem, alternate=None):
    # Regression Problem
    if problem == "regression":
        return Regression(alternate)
    elif problem == "classification":
        return Regression(alternate, error=inverse_f1_score)
    # String Match Problem
    elif problem == "string_match":
        return StringMatch(alternate)
    elif problem == "progsys":
        return ProgSys(alternate)
    elif problem == "new":
        print("new problem goes here")
        # parameters.FITNESS_FUNCTION = whatever
    else:
        print("Please specify a valid fitness function")
        exit()
