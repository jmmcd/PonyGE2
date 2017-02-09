from datetime import datetime
from time import time
from random import seed
from sys import version_info
from socket import gethostname
from os import getpid

from algorithm.parameters import params
from stats.stats import generate_folders_and_files
from utilities.stats import trackers


def check_python_version():
    """
    Check the python version to ensure it is correct. PonyGE uses Python 3.

    :return: Nothing
    """

    if version_info.major < 3 or version_info.minor < 5:
        print("\nError: Python version not supported. Must use at least "
              "Python 3.5")
        quit()


def initialise_run_params():
    """
    Initialises all lists and trackers. Generates save folders and initial
    parameter files if debugging is not active.

    :return: Nothing
    """

    start = datetime.now()
    trackers.time_list.append(time())

    # Set random seed
    if params['RANDOM_SEED'] is None:
        params['RANDOM_SEED'] = int(start.microsecond)
    seed(params['RANDOM_SEED'])

    # Generate a time stamp for use with folder and file names.
    hms = "%02d%02d%02d" % (start.hour, start.minute, start.second)
    params['TIME_STAMP'] = "_".join([gethostname(),
                                     str(params['RANDOM_SEED']),
                                     str(getpid()),
                                     str(start.year)[2:],
                                     str(start.month),
                                     str(start.day),
                                     hms,
                                     str(start.microsecond)])
    if not params['SILENT']:
        print("\nStart:\t", start, "\n")

    # Generate save folders and files
    if params['DEBUG']:
        print("Seed:\t", params['RANDOM_SEED'], "\n")
    else:
        generate_folders_and_files()


def make_import_str(fns, location):
    """
    Takes in a paired list of operators and the specified function from the
    option parser. Strings either represent the full dotted path to the
    function which we wish to access, eg operators.selection.tournament, or
    just the function name directly (in which case we default to the specified
    functions in the default location for each operators).

    :param fns: a paired list of operators and the specified function from the
    option parser. Strings either represent the full dotted path to the
    function which we wish to access, eg operators.selection.tournament, or
    just the function name directly (in which case we default to the specified
    functions in the default location for each operators).
    :param location: A string specifying the containing folder of the
    functions listed in fns, e.g. "operators", "fitness", "utilities", etc.
    :return: a string of imports of correct modules,
    eg import operators.selection
    """

    imports = []
    for fn in [func for func in fns if func[1]]:
        # Extract pairs of operators and functions, but only if functions
        # are not 'None' (some default parameters e.g. error_metric may be
        # set to 'None').
        operator, function = fn[0], fn[1]
        parts = function.split(".")
        # Split the function into its component parts

        if len(parts) == 1:
            # If the specified location is a single name, default to
            # operators.operator location
            imports.append("import " + location + "." + operator.lower())
            params[operator] = ".".join([location, operator.lower(), parts[0]])

        else:
            # "operators.selection.tournament" -> "import operators.selection"
            imports.append("import " + ".".join(parts[:-1]))

    return "\n".join(imports)


def set_param_imports():
    """
    This function makes the command line experience easier for users. When
    specifying operators listed in the lists below, users do not need to
    specify the full file path to the functions themselves. Users can simply
    specify a single word, e.g.

        "--mutation subtree"

    Using the special_ops dictionary for example, this will default to
    "operators.mutation.subtree. Executes the correct imports for specified
    modules and then saves the correct parameters in the params dictionary.
    Users can still specify the full direct path to the operators if they so
    desire, allowing them to create new operators and save them wherever
    they like.

    Sets the fitness function for a problem automatically. Fitness functions
    are stored in fitness. Fitness functions must be classes, where the
    class name matches the file name.

    Function is set up to automatically set imports for operators and error
    metrics.

    :return: Nothing.
    """

    # For these ops we let the param equal the function itself.
    ops = {'operators': ['INITIALISATION', 'SELECTION', 'CROSSOVER',
                         'MUTATION', 'REPLACEMENT'],
           'utilities.fitness': ['ERROR_METRIC'],
           'fitness': ['FITNESS_FUNCTION'],
           'algorithm': ['SEARCH_LOOP', 'STEP']}

    # We have to take 'algorithm' first as the functions from
    # algorithm need to be imported before any others to prevent
    # circular imports. We have to take 'utilities.fitness' before
    # 'fitness' because ERROR_METRIC has to be set in order to call
    # the fitness function constructor.

    for special_ops in ['algorithm', 'utilities.fitness', 'operators', 'fitness']:

        if all([callable(params[op]) for op in ops[special_ops]]):
            # params are already functions
            pass

        else:
            if special_ops == "fitness":
                import_func = "from fitness." + params[
                    'FITNESS_FUNCTION'] + " import " + params[
                                  'FITNESS_FUNCTION']

                # Import the required fitness function.
                exec(import_func)

                # Set the fitness function in the params dictionary.
                params['FITNESS_FUNCTION'] = eval(params['FITNESS_FUNCTION'] +
                                                  "()")
            else:
                # We need to do an appropriate import...
                import_str = make_import_str([[op, params[op]] for op in
                                              ops[special_ops]], special_ops)

                exec(import_str)
                # ... and then eval the param.
                for op in ops[special_ops]:
                    if params[op]:
                        params[op] = eval(params[op])
