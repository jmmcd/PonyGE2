import argparse
import importlib
import numpy as np

from datetime import datetime
from math import ceil
from operator import attrgetter
from os import getpid, path
from platform import system

from parameters import params
from representation import Grammar
from operators import load_population
from random import seed
from datetime import datetime
from socket import gethostname
from time import time
from os import getpid, path
import trackers
from arc_evaluate import arc_evaluate

"""
command_line_parser.py
"""

class SortingHelpFormatter(argparse.HelpFormatter):
    """
    Custom class for sorting the arguments of the arg parser for printing. When
    "--help" is called, arguments will be listed in alphabetical order. Without
    this custom class, arguments will be printed in the order in which they are
    defined.
    """
    
    def add_arguments(self, actions):
        actions = sorted(actions, key=attrgetter('option_strings'))
        super(SortingHelpFormatter, self).add_arguments(actions)

def parse_cmd_args(arguments):
    """
    Parser for command line arguments specified by the user. Specified command
    line arguments over-write parameter file arguments, which themselves
    over-write original values in the algorithm.parameters.params dictionary.

    The argument parser structure is set up such that each argument has the
    following information:

        dest: a valid key from the algorithm.parameters.params dictionary
        type: an expected type for the specified option (i.e. str, int, float)
        help: a string detailing correct usage of the parameter in question.

    Optional info:

        default: The default setting for this parameter.
        action : The action to be undertaken when this argument is called.

    NOTE: You cannot add a new parser argument and have it evaluate "None" for
    its value. All parser arguments are set to "None" by default. We filter
    out arguments specified at the command line by removing any "None"
    arguments. Therefore, if you specify an argument as "None" from the
    command line and you evaluate the "None" string to a None instance, then it
    will not be included in the eventual parameters.params dictionary. A
    workaround for this would be to leave "None" command line arguments as
    strings and to eval them at a later stage.

    :param arguments: Command line arguments specified by the user.
    :return: A dictionary of parsed command line arguments, along with a
    dictionary of newly specified command line arguments which do not exist
    in the params dictionary.
    """

    # Initialise parser
    parser = argparse.ArgumentParser(
        formatter_class=SortingHelpFormatter,
        usage=argparse.SUPPRESS,
        description="""Welcome to PonyGE2 - Help.
        The following are the available command line arguments. Please see
        src/algorithm/parameters.py for a more detailed explanation of each
        argument and its possible values.""",
        epilog="""To try out PonyGE2 from the command line simply navigate to
        the src directory and type: python ponyge.py.""")

    parser._optionals.title = 'PonyGE2 command-line usage'

    # Set up class for parsing list arguments.
    class ListAction(argparse.Action):
        """
        Class for parsing a given string into a list.
        """

        def __init__(self, option_strings, **kwargs):
            super(ListAction, self).__init__(option_strings, **kwargs)
    
        def __call__(self, parser, namespace, value, option_string=None):
            if type(eval(value)) != list or any([type(i) != int for i in
                                                 eval(value)]):
                s = "utilities.algorithm.command_line_parser.ListAction\n" \
                    "Error: parameter %s is not a valid genome.\n" \
                    "       Value given: %s" % (option_string, value)
                raise Exception(s)
            else:
                setattr(namespace, self.dest, eval(value))

    # Set up class for checking float arguments.
    class FloatAction(argparse.Action):
        """
        Class for checking a given float is within the range [0:1].
        """

        def __init__(self, option_strings, **kwargs):
            super(FloatAction, self).__init__(option_strings, **kwargs)

        def __call__(self, parser, namespace, value, option_string=None):
            if not 0 <= float(value) <= 1:
                s = "utilities.algorithm.command_line_parser.FloatAction\n" \
                    "Error: parameter %s outside allowed range [0:1].\n" \
                    "       Value given: %s" % (option_string, value)
                raise Exception(s)
            else:
                setattr(namespace, self.dest, float(value))

    # Set up class for checking raw string arguments to catch "tab" inputs.
    class CatchTabStr(argparse.Action):
        """
        Class for checking raw string arguments to catch "tab" inputs.
        """
    
        def __init__(self, option_strings, **kwargs):
            super(CatchTabStr, self).__init__(option_strings, **kwargs)
    
        def __call__(self, parser, namespace, value, option_string=None):
            if repr(value) == repr("\\t"):
                value = "\t"
            setattr(namespace, self.dest, value)

    # LOAD PARAMETERS FILE
    parser.add_argument('--parameters',
                        dest='PARAMETERS',
                        type=str,
                        help='Specifies the parameters file to be used. Must '
                             'include the full file extension. Full file path'
                             'does NOT need to be specified.')

    # LOAD STEP AND SEARCH LOOP FUNCTIONS
    parser.add_argument('--search_loop',
                        dest='SEARCH_LOOP',
                        type=str,
                        help='Sets the desired search loop function.')
    parser.add_argument('--step',
                        dest='STEP',
                        type=str,
                        help='Sets the desired search step function.')

    # POPULATION OPTIONS
    parser.add_argument('--population_size',
                        dest='POPULATION_SIZE',
                        type=int,
                        help='Sets the population size, requires int value.')
    parser.add_argument('--generations',
                        dest='GENERATIONS',
                        type=int,
                        help='Sets the number of generations, requires int '
                             'value.')
    parser.add_argument('--hill_climbing_history',
                        dest='HILL_CLIMBING_HISTORY',
                        type=int,
                        help='Sets the history-length for late-acceptance'
                        'and step-counting hill-climbing.')
    parser.add_argument('--schc_count_method',
                        dest='SCHC_COUNT_METHOD',
                        type=str,
                        help='Sets the counting method for step-counting '
                             'hill-climbing. Optional values are "count_all", '
                             '"acp", and "imp".')

    # INDIVIDUAL SIZE
    parser.add_argument('--max_tree_depth',
                        dest='MAX_TREE_DEPTH',
                        type=int,
                        help='Sets the max derivation tree depth for the '
                             'algorithm, requires int value. The default max '
                             'tree depth is set to None, i.e. trees can grow'
                             'indefinitely. This can also be set by '
                             'specifying the max tree depth to be 0.')
    parser.add_argument('--max_tree_nodes',
                        dest='MAX_TREE_NODES',
                        type=int,
                        help='Sets the max derivation tree nodes for the '
                             'algorithm, requires int value. The default max '
                             'tree nodes is set to None, i.e. trees can grow'
                             'indefinitely. This can also be set by '
                             'specifying the max tree nodes to be 0.')
    parser.add_argument('--codon_size',
                        dest='CODON_SIZE',
                        type=int,
                        help='Sets the range from 0 to condon_size to be used '
                             'in genome, requires int value')
    parser.add_argument('--max_genome_length',
                        dest='MAX_GENOME_LENGTH',
                        type=int,
                        help='Sets the maximum chromosome length for the '
                             'algorithm, requires int value. The default max '
                             'genome length is set to None, i.e. gemomes can '
                             'grow indefinitely. This can also be set by '
                             'specifying the max genome length to be 0.')
    parser.add_argument('--max_wraps',
                        dest='MAX_WRAPS',
                        type=int,
                        help='Sets the maximum number of times the genome '
                             'mapping process can wrap over the length of the '
                             'genome. Requires int value.')
    parser.add_argument('--permutation_ramps',
                        dest='PERMUTATION_RAMPS',
                        type=int,
                        help='Set the number of depths permutations are '
                             'calculated for (starting from the minimum path '
                             'of the grammar). Mainly for use with '
                             'the grammar analyser script. Requires int '
                             'value.')

    # INITIALISATION
    parser.add_argument('--max_init_tree_depth',
                        dest='MAX_INIT_TREE_DEPTH',
                        type=int,
                        help='Sets the max tree depth for initialisation.')
    parser.add_argument('--min_init_tree_depth',
                        dest='MIN_INIT_TREE_DEPTH',
                        type=int,
                        help='Sets the min tree depth for initialisation.')
    parser.add_argument('--init_genome_length',
                        dest='INIT_GENOME_LENGTH',
                        type=int,
                        help='Sets the length for chromosomes to be '
                             'initialised to. Requires int value.')
    parser.add_argument('--initialisation',
                        dest='INITIALISATION',
                        type=str,
                        help='Sets the initialisation strategy, requires a '
                             'string such as "rhh" or a direct path string '
                             'such as "operators.initialisation.rhh".')

    # SELECTION
    parser.add_argument('--selection',
                        dest='SELECTION',
                        type=str,
                        help='Sets the selection to be used, requires string '
                             'such as "tournament" or direct path string such '
                             'as "operators.selection.tournament".')
    parser.add_argument('--invalid_selection',
                        dest='INVALID_SELECTION',
                        action='store_true',
                        default=None,
                        help='Allow for the selection of invalid individuals '
                             'during selection.')
    parser.add_argument('--tournament_size',
                        dest='TOURNAMENT_SIZE',
                        type=int,
                        help='Sets the number of indivs to contest tournament,'
                             ' requires int.')
    parser.add_argument('--selection_proportion',
                        dest='SELECTION_PROPORTION',
                        action=FloatAction,
                        help='Sets the proportion for truncation selection, '
                             'requires float, e.g. 0.5.')

    # OPERATOR OPTIONS
    parser.add_argument('--within_used',
                        dest='WITHIN_USED',
                        default=None,
                        action='store_true',
                        help='Boolean flag for selecting whether or not '
                             'mutation is confined to within the used portion '
                             'of the genome. Default set to True.')
    
    # CROSSOVER
    parser.add_argument('--crossover',
                        dest='CROSSOVER',
                        type=str,
                        help='Sets the type of crossover to be used, requires '
                             'string such as "subtree" or direct path string '
                             'such as "operators.crossover.subtree".')
    parser.add_argument('--crossover_probability',
                        dest='CROSSOVER_PROBABILITY',
                        action=FloatAction,
                        help='Sets the crossover probability, requires float, '
                             'e.g. 0.9.')
    parser.add_argument('--no_crossover_invalids',
                        dest='NO_CROSSOVER_INVALIDS',
                        default=None,
                        action='store_true',
                        help='Prevents invalid individuals from being '
                             'generated by crossover.')

    # MUTATION
    parser.add_argument('--mutation',
                        dest='MUTATION',
                        type=str,
                        help='Sets the type of mutation to be used, requires '
                             'string such as "int_flip_per_codon" or direct '
                             'path string such as '
                             '"operators.mutation.int_flip_per_codon".')
    parser.add_argument('--mutation_events',
                        dest='MUTATION_EVENTS',
                        type=int,
                        help='Sets the number of mutation events based on '
                             'probability.')
    parser.add_argument('--mutation_probability',
                        dest='MUTATION_PROBABILITY',
                        action=FloatAction,
                        help='Sets the rate of mutation probability for linear'
                             ' genomes')
    parser.add_argument('--no_mutation_invalids',
                        dest='NO_MUTATION_INVALIDS',
                        default=None,
                        action='store_true',
                        help='Prevents invalid individuals from being '
                             'generated by mutation.')

    # EVALUATION
    parser.add_argument('--fitness_function',
                        dest='FITNESS_FUNCTION',
                        type=str,
                        nargs='+',
                        help='Sets the fitness function to be used. '
                             'Requires string such as "regression". '
                             'Multiple fitness functions can be specified'
                             'for multiple objective optimisation (using '
                             'NSGA-II). To specify multiple fitness '
                             'functions simply enter in the desired names of'
                             ' the functions separated by spaces.')
    parser.add_argument('--dataset_train',
                        dest='DATASET_TRAIN',
                        type=str,
                        help='For use with problems that use a dataset. '
                             'Specifies the training data for evolution. '
                             'Full file name must be specified.')
    parser.add_argument('--dataset_test',
                        dest='DATASET_TEST',
                        type=str,
                        help='For use with problems that use a dataset. '
                             'Specifies the testing data for evolution. '
                             'Full file name must be specified.')
    parser.add_argument('--dataset_delimiter',
                        dest='DATASET_DELIMITER',
                        action=CatchTabStr,
                        help='For use with problems that use a dataset. '
                             'Specifies the delimiter for the dataset. '
                             'Requires string such as "\\t".')
    parser.add_argument('--target',
                        dest='TARGET',
                        type=str,
                        help='For string match problem. Requires target '
                             'string.')
    parser.add_argument('--error_metric',
                        dest='ERROR_METRIC',
                        type=str,
                        help='Sets the error metric to be used with supervised'
                             ' learning problems. Requires string such as '
                             '"mse" or "rmse".')
    parser.add_argument('--optimize_constants',
                        dest='OPTIMIZE_CONSTANTS',
                        action='store_true',
                        default=None,
                        help='Whether to optimize numerical constants by '
                             'gradient descent in supervised learning '
                             'problems. Requires True or False, default '
                             'False.')
    parser.add_argument('--multicore',
                        dest='MULTICORE',
                        action='store_true',
                        default=None,
                        help='Turns on multicore evaluation.')
    parser.add_argument('--cores',
                        dest='CORES',
                        type=int,
                        help='Specify the number of cores to be used for '
                             'multicore evaluation. Requires int.')
    
    # REPLACEMENT
    parser.add_argument('--replacement',
                        dest='REPLACEMENT',
                        type=str,
                        help='Sets the replacement strategy, requires string '
                             'such as "generational" or direct path string '
                             'such as "operators.replacement.generational".')
    parser.add_argument('--elite_size',
                        dest='ELITE_SIZE',
                        type=int,
                        help='Sets the number of elites to be used, requires '
                             'int value.')

    # PROBLEM SPECIFICS
    parser.add_argument('--grammar_file',
                        dest='GRAMMAR_FILE',
                        type=str,
                        help='Sets the grammar to be used, requires string.')
    parser.add_argument('--experiment_name',
                        dest='EXPERIMENT_NAME',
                        type=str,
                        help='Optional parameter to save results in '
                             'results/[EXPERIMENT_NAME] folder. If not '
                             'specified then results are saved in default '
                             'results folder.')
    parser.add_argument('--runs',
                        dest='RUNS',
                        type=int,
                        help='Optional parameter to specify the number of '
                             'runs to be performed for an experiment. Only '
                             'used with experiment manager.')
    parser.add_argument('--extra_parameters',
                        dest='EXTRA_PARAMETERS',
                        type=str,
                        nargs='+',
                        help='Optional extra command line parameter for '
                             'inclusion of any extra information required '
                             'for user-specific runs. Can be whatever you '
                             'want it to be. Specified arguments are parsed '
                             'as a list. Specify as many values as desired, '
                             'separated by spaces.')

    # OPTIONS
    parser.add_argument('--random_seed',
                        dest='RANDOM_SEED',
                        type=int,
                        help='Sets the random seed to be used with both the '
                             'standard Python RNG and the NumPy RNG. '
                             'requires int value.')
    parser.add_argument('--debug',
                        dest='DEBUG',
                        action='store_true',
                        default=None,
                        help='Disables saving of all ancillary files.')
    parser.add_argument('--verbose',
                        dest='VERBOSE',
                        action='store_true',
                        default=None,
                        help='Turns on the verbose output of the program in '
                             'terms of command line and extra files.')
    parser.add_argument('--silent',
                        dest='SILENT',
                        action='store_true',
                        default=None,
                        help='Prevents any output from being printed to the '
                             'command line.')
    parser.add_argument('--save_all',
                        dest='SAVE_ALL',
                        action='store_true',
                        default=None,
                        help='Saves the best phenotypes at each generation.')
    parser.add_argument('--save_plots',
                        dest='SAVE_PLOTS',
                        action='store_true',
                        default=None,
                        help='Saves plots for best fitness.')

    # REVERSE-MAPPING
    parser.add_argument('--reverse_mapping_target',
                        dest='REVERSE_MAPPING_TARGET',
                        type=str,
                        help='Target string to parse into a GE individual.')
    parser.add_argument('--target_seed_folder',
                        dest='TARGET_SEED_FOLDER',
                        type=str,
                        help='Specify a target seed folder in the "seeds" '
                             'directory that contains a population of '
                             'individuals with which to seed a run.')

    # STATE SAVING/LOADING
    parser.add_argument('--save_state',
                        dest='SAVE_STATE',
                        action='store_true',
                        default=None,
                        help='Saves the state of the evolutionary run every '
                             'generation. You can specify how often you want '
                             'to save the state with the command '
                             '"--save_state_step".')
    parser.add_argument('--save_state_step',
                        dest='SAVE_STATE_STEP',
                        type=int,
                        help='Specifies how often the state of the current '
                             'evolutionary run is saved (i.e. every n-th '
                             'generation). Requires int value.')
    parser.add_argument('--load_state',
                        dest='LOAD_STATE',
                        type=str,
                        help='Load an evolutionary run from a saved state. '
                             'You must specify the full file path to the '
                             'desired state file. Note that state files have '
                             'no file type.')


    # MULTIAGENT
    parser.add_argument('--multiagent',
                        dest='MULTIAGENT',
                        action='store_true',
                        default=None,                        
                        help='This enable the multiagent mode. If this mode is'
                             ' enabled the search_loop and step parameter are' 
                             ' overridden with search_multiagent and step_multiagent'
                             ' respectively')
    parser.add_argument('--agent_size',
                        dest='AGENT_SIZE',
                        type=int,
                        help='Specifies how many agents are initialize in'
                             ' the environment. By default 100 agents are initialize.'
                             ' Greater the number of agents the time to find the'
                             ' would be reduced')
    parser.add_argument('--interaction_probability',
                        dest='INTERACTION_PROBABILITY',
                        action=FloatAction,
                        help='Specifies the probability of agent interacting with'
                             ' other nearby agents in the environment. By default'
                             ' 0.5 probability is used. Higher the probability the time'
                             ' to find the solution would be reduced')

    # CACHING
    class CachingAction(argparse.Action):
        """
        Class for defining special mutually exclusive options for caching.
        """

        def __init__(self, option_strings, CACHE=None, LOOKUP_FITNESS=None,
                     LOOKUP_BAD_FITNESS=None, MUTATE_DUPLICATES=None,
                     **kwargs):
            self.CACHE = CACHE
            self.LOOKUP_FITNESS = LOOKUP_FITNESS
            self.LOOKUP_BAD_FITNESS = LOOKUP_BAD_FITNESS
            self.MUTATE_DUPLICATES = MUTATE_DUPLICATES
            super(CachingAction, self).__init__(option_strings, nargs=0,
                                                **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            setattr(namespace, 'CACHE', self.CACHE)
            if 'LOOKUP_FITNESS' not in namespace or \
                            getattr(namespace, 'LOOKUP_FITNESS') is not False:
                # able to overwrite if True or None
                setattr(namespace, 'LOOKUP_FITNESS', self.LOOKUP_FITNESS)
            if self.LOOKUP_BAD_FITNESS and \
                            'LOOKUP_BAD_FITNESS' not in namespace:
                setattr(namespace, 'LOOKUP_BAD_FITNESS',
                        self.LOOKUP_BAD_FITNESS)
            if self.MUTATE_DUPLICATES and 'MUTATE_DUPLICATES' not in namespace:
                setattr(namespace, 'MUTATE_DUPLICATES', self.MUTATE_DUPLICATES)

    # Generate a mutually exclusive group for caching options. This means
    # that you cannot specify multiple caching options simultaneously,
    # only one at a time.
    parser.add_argument("--cache",
                        dest='CACHE',
                        action=CachingAction,
                        CACHE=True,
                        LOOKUP_FITNESS=True,
                        help='Tracks unique phenotypes and is used to '
                             'lookup duplicate fitnesses.')
    caching_group = parser.add_mutually_exclusive_group()
    caching_group.add_argument("--dont_lookup_fitness",
                               dest='CACHE',
                               action=CachingAction,
                               CACHE=True,
                               LOOKUP_FITNESS=False,
                               help='Uses cache to track duplicate '
                                    'individuals, but does not use the cache '
                                    'to save fitness evaluations.')
    caching_group.add_argument("--lookup_bad_fitness",
                               dest='CACHE',
                               action=CachingAction,
                               CACHE=True,
                               LOOKUP_FITNESS=False,
                               LOOKUP_BAD_FITNESS=True,
                               help='Gives duplicate phenotypes a bad fitness '
                                    'when encountered. Uses cache.')
    caching_group.add_argument("--mutate_duplicates",
                               dest='CACHE',
                               action=CachingAction,
                               CACHE=True,
                               LOOKUP_FITNESS=False,
                               MUTATE_DUPLICATES=True,
                               help='Replaces duplicate individuals with '
                                    'mutated versions. Uses cache.')

    # Parse command line arguments using all above information.
    args, unknown = parser.parse_known_args(arguments)

    # All default args in the parser are set to "None". Only take arguments
    # which are not "None", i.e. arguments which have been passed in from
    # the command line.
    cmd_args = {key: value for key, value in vars(args).items() if value is
                not None}

    # Set "None" values correctly.
    for key in sorted(cmd_args.keys()):
        # Check all specified arguments.
        
        if type(cmd_args[key]) == str and cmd_args[key].lower() == "none":
            # Allow for people not using correct capitalisation.
            
            cmd_args[key] = None

    return cmd_args, unknown

"""
math_functions.py
"""

def return_one_percent(num, pop_size):
    """
    Returns either one percent of the population size or a given number,
    whichever is larger.

    :param num: A given number of individuals (NOT a desired percentage of
    the population).
    :param pop_size: A given population size.
    :return: either one percent of the population size or a given number,
    whichever is larger.
    """

    # Calculate one percent of the given population size.
    percent = int(round(pop_size/100))

    # Return the biggest number.
    if percent < num:
        return num
    else:
        return percent

"""
initialise_run.py
"""

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
    must be stored in fitness. Fitness functions must be classes, where the
    class name matches the file name.

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

    for special_ops in ['algorithm', 'utilities.fitness',
                        'operators', 'fitness']:

        if all([callable(params[op]) for op in ops[special_ops]]):
            # params are already functions
            pass

        else:

            for op in ops[special_ops]:

                if special_ops == "fitness":
                    # Fitness functions represent a special case.
                    get_fit_func_imports()

                elif params[op] is not None:
                    # Split import name based on "." to find nested modules.
                    split_name = params[op].split(".")

                    if len(split_name) > 1:
                        # Check to see if a full path has been specified.

                        # Get attribute name.
                        attr_name = split_name[-1]

                        try:
                            # Try and use the exact specified path to load
                            # the module.

                            # Get module name.
                            module_name = ".".join(split_name[:-1])

                            # Import module and attribute and save.
                            params[op] = return_attr_from_module(module_name,
                                                                 attr_name)

                        except Exception:
                            # Either a full path has not actually been
                            # specified, or the module doesn't exist. Try to
                            # append specified module to default location.

                            # Get module name.
                            module_name = ".".join([special_ops,
                                                    ".".join(split_name[:-1])])

                            try:
                                # Import module and attribute and save.
                                params[op] = return_attr_from_module(module_name,
                                                                     attr_name)

                            except Exception:
                                s = "utilities.algorithm.initialise_run." \
                                    "set_param_imports\n" \
                                    "Error: Specified %s function not found:" \
                                    " %s\n" \
                                    "       Checked locations: %s\n" \
                                    "                          %s\n" \
                                    "       Please ensure parameter is " \
                                    "specified correctly." % \
                                    (op.lower(), attr_name, params[op],
                                     ".".join([module_name, attr_name]))
                                raise Exception(s)

                    else:
                        # Just module name specified. Use default location.

                        # If multiagent is specified need to change
                        # how search and step module is called
                        # Loop and step functions for multiagent is contained 
                        # inside algorithm search_loop_distributed and 
                        # step_distributed respectively

                        if params['MULTIAGENT'] and \
                        ( op == 'SEARCH_LOOP' or op == 'STEP' ) :
                            # Define the directory structure for the multiagent search
                            # loop and step
                            multiagent_ops = {'search_loop':'distributed_algorithm.search_loop' \
                                                ,'step':'distributed_algorithm.step'}

                            # Get module and attribute names
                            module_name = ".".join([special_ops, multiagent_ops[op.lower()]])
                            attr_name = split_name[-1]

                        else:
                            # Get module and attribute names.
                            module_name = ".".join([special_ops, op.lower()])
                            attr_name = split_name[-1]

                        # Import module and attribute and save.
                        params[op] = return_attr_from_module(module_name,
                                                             attr_name)

def get_fit_func_imports():
    """
    Special handling needs to be done for fitness function imports,
    as fitness functions can be specified a number of different ways. Notably,
    a list of fitness functions can be specified, indicating multiple
    objective optimisation.

    Note that fitness functions must be classes where the class has the same
    name as its containing file. Fitness functions must be contained in the
    `fitness` module.

    :return: Nothing.
    """

    op = 'FITNESS_FUNCTION'

    if "," in params[op]:
        # List of fitness functions given in parameters file.

        # Convert specified fitness functions into a list of strings.
        params[op] = params[op].strip("[()]").split(",")

    if isinstance(params[op], list) and len(params[op]) == 1:
        # Single fitness function given in a list format. Don't use
        # multi-objective optimisation.

        params[op] = params[op][0]

    if isinstance(params[op], list):
        # List of multiple fitness functions given.

        for i, name in enumerate(params[op]):

            # Split import name based on "." to find nested modules.
            split_name = name.strip().split(".")

            # Get module and attribute names.
            module_path = ".".join(['fitness', name.strip()])
            attr = split_name[-1]

            # Import this fitness function.
            params[op][i] = return_attr_from_module(module_path, attr)

        # # Import base multi-objective fitness function class.
        # from fitness.base_ff_classes.moo_ff import moo_ff

        # # Set main fitness function as base multi-objective fitness
        # # function class.
        # params[op] = moo_ff(params[op])

    else:
        # A single fitness function has been specified.
        if params[op] == "arc_evaluate":
            params[op] = arc_evaluate()
        else:
            # Split import name based on "." to find nested modules.
            split_name = params[op].strip().split(".")

            # Get attribute name.
            attr_name = split_name[-1]

            # Get module name.
            module_name = ".".join(["fitness", params[op]])

            # Import module and attribute and save.
            params[op] = return_attr_from_module(module_name, attr_name)

            # Initialise fitness function.
            params[op] = params[op]()

def return_attr_from_module(module_name, attr_name):
    """
    Given a module path and the name of an attribute that exists in that
    module, import the attribute from the module using the importlib package
    and return it.

    :param module_name: The name/location of the desired module.
    :param attr_name: The name of the attribute.
    :return: The imported attribute from the module.
    """

    try:
        # Import module.
        module = importlib.import_module(module_name)

    except ModuleNotFoundError:
        s = "utilities.algorithm.initialise_run.return_attr_from_module\n" \
            "Error: Specified module not found: %s" % (module_name)
        raise Exception(s)

    try:
        # Import specified attribute and return.
        return getattr(module, attr_name)

    except AttributeError:
        s = "utilities.algorithm.initialise_run.return_attr_from_module\n" \
            "Error: Specified attribute '%s' not found in module '%s'." \
            % (attr_name, module_name)
        raise Exception(s)

def pool_init(params_):
    """
    When initialising the pool the original params dict (params_) is passed in
    and used to update the newly created instance of params, as Windows does
    not retain the system memory of the parent process.

    :param params_: original params dict
    :return: Nothing.
    """


    if system() == 'Windows':
        params.update(params_)

def set_params(command_line_args, create_files=True):
    """
    This function parses all command line arguments specified by the user.
    If certain parameters are not set then defaults are used (e.g. random
    seeds, elite size). Sets the correct imports given command line
    arguments. Sets correct grammar file and fitness function. Also
    initialises save folders and tracker lists in utilities.trackers.

    :param command_line_args: Command line arguments specified by the user.
    :return: Nothing.
    """
    cmd_args, unknown = parse_cmd_args(command_line_args)

    if unknown:
        # We currently do not parse unknown parameters. Raise error.
        s = "algorithm.parameters.set_params\nError: " \
            "unknown parameters: %s\nYou may wish to check the spelling, " \
            "add code to recognise this parameter, or use " \
            "--extra_parameters" % str(unknown)
        raise Exception(s)


    # Join original params dictionary with command line specified arguments.
    # NOTE that command line arguments overwrite all previously set parameters.
    params.update(cmd_args)

    if params['REPLACEMENT'].split(".")[-1] == "steady_state":
        # Set steady state step and replacement.
        params['STEP'] = "steady_state_step"
        params['GENERATION_SIZE'] = 2

    else:
        # Elite size is set to either 1 or 1% of the population size,
        # whichever is bigger if no elite size is previously set.
        if params['ELITE_SIZE'] is None:
            params['ELITE_SIZE'] = return_one_percent(1, params[
                'POPULATION_SIZE'])

        # Set the size of a generation
        params['GENERATION_SIZE'] = params['POPULATION_SIZE'] - \
                                    params['ELITE_SIZE']

    # Initialise run lists and folders before we set imports.r
    initialise_run_params(create_files)

    # Set correct param imports for specified function options, including
    # error metrics and fitness functions.
    set_param_imports()

    # Set GENOME_OPERATIONS automatically for faster linear operations.
    if (params['CROSSOVER'].representation == "subtree" or
        params['MUTATION'].representation == "subtree"):
        params['GENOME_OPERATIONS'] = False
    else:
        params['GENOME_OPERATIONS'] = True

    # Ensure correct operators are used if multiple fitness functions used.
    if hasattr(params['FITNESS_FUNCTION'], 'multi_objective'):

        # Check that multi-objective compatible selection is specified.
        if not hasattr(params['SELECTION'], "multi_objective"):
            s = "algorithm.parameters.set_params\n" \
                "Error: multi-objective compatible selection " \
                "operator not specified for use with multiple " \
                "fitness functions."
            raise Exception(s)

        if not hasattr(params['REPLACEMENT'], "multi_objective"):

            # Check that multi-objective compatible replacement is
            # specified.
            if not hasattr(params['REPLACEMENT'], "multi_objective"):
                s = "algorithm.parameters.set_params\n" \
                    "Error: multi-objective compatible replacement " \
                    "operator not specified for use with multiple " \
                    "fitness functions."
                raise Exception(s)

    # Parse grammar file and set grammar class.
    params['BNF_GRAMMAR'] = Grammar(path.join("..", "grammars",
                                            params['GRAMMAR_FILE']))

    # Population loading for seeding runs (if specified)
    if params['TARGET_SEED_FOLDER']:


        # A target folder containing seed individuals has been given.
        params['SEED_INDIVIDUALS'] = load_population(
            params['TARGET_SEED_FOLDER'])

def initialise_run_params(create_files):
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
                                     str(start.year)[2:],
                                     str(start.month),
                                     str(start.day),
                                     hms,
                                     str(start.microsecond),
                                     str(getpid()),
                                     str(params['RANDOM_SEED'])])
