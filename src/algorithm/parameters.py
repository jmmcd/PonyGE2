from multiprocessing import cpu_count
from os import path
from socket import gethostname

hostname = gethostname().split('.')
machine_name = hostname[0]


"""Algorithm parameters"""
params = {
        # Set default step and search loop functions
        'SEARCH_LOOP': 'search_loop',
        'STEP': 'step',

        # Evolutionary Parameters
        'POPULATION_SIZE': 500,
        'GENERATIONS': 50,
        'HILL_CLIMBING_HISTORY': 1000,

        # Set optional experiment name
        'EXPERIMENT_NAME': None,

        # Class of problem
        'FITNESS_FUNCTION': "regression",
        # "regression"
        # "string_match"
        # "classification"

        # Select problem dataset
        'DATASET': "Vladislavleva4",
        # "Dow"
        # "Keijzer6"
        # "Vladislavleva4"

        # Set grammar file
        'GRAMMAR_FILE': "Vladislavleva4.bnf",
        # "Vladislavleva4.bnf"
        # "Keijzer6.bnf"
        # "Dow.bnf"
        # "Banknote.bnf"
        # "letter.bnf"

        # Select error metric
        'ERROR_METRIC': None,
        # "mse"
        # "mae"
        # "rmse"
        # "hinge"
        # "f1_score"

        # Specify target for target problems
        'TARGET': "ponyge_rocks",

        # Set max sizes of individuals
        'MAX_TREE_DEPTH': 17,
        'MAX_TREE_NODES': None,
        'CODON_SIZE': 100000,
        'MAX_INIT_GENOME_LENGTH': 200,
        'MAX_GENOME_LENGTH': 500,
        'MAX_WRAPS': 0,

        # INITIALISATION
        'INITIALISATION': "operators.initialisation.rhh",
        # "operators.initialisation.uniform_genome"
        # "operators.initialisation.rhh"
        'MAX_INIT_TREE_DEPTH': 10,
        # Set the maximum tree depth for initialisation.
        'MIN_INIT_TREE_DEPTH': None,
        # Set the minimum tree depth for initialisation.

        # SELECTION
        'SELECTION': "operators.selection.tournament",
        # "operators.selection.tournament"
        # "operators.selection.truncation",
        'TOURNAMENT_SIZE': 2,
        # For tournament selection
        'SELECTION_PROPORTION': 0.5,
        # For truncation selection
        'INVALID_SELECTION': False,
        # Allow for selection of invalid individuals during selection process.

        # CROSSOVER
        'CROSSOVER': "operators.crossover.variable_onepoint",
        # "operators.crossover.fixed_onepoint",
        # "operators.crossover.subtree",
        'CROSSOVER_PROBABILITY': 0.75,
        'NO_CROSSOVER_INVALIDS': False,
        # Prevents crossover from generating invalids.

        # MUTATION
        'MUTATION': "operators.mutation.int_flip",
        # "operators.mutation.subtree",
        # "operators.mutation.int_flip",
        'MUTATION_PROBABILITY': None,
        'MUTATION_EVENTS': 1,
        'NO_MUTATION_INVALIDS': False,
        # Prevents mutation from generating invalids.

        # REPLACEMENT
        'REPLACEMENT': "operators.replacement.generational",
        # "operators.replacement.generational",
        # "operators.replacement.steady_state",
        'ELITE_SIZE': None,

        # DEBUGGING
        # Use this to turn on debugging mode. This mode doesn't write any files
        # and should be used when you want to test new methods.
        'DEBUG': False,

        # PRINTING
        # Use this to print out basic statistics for each generation to the
        # command line.
        'VERBOSE': False,
        # Use this to prevent anything being printed to the command line.
        'SILENT': False,

        # SAVING
        'SAVE_ALL': False,
        # Use this to save the phenotype of the best individual from each
        # generation. Can generate a lot of files. DEBUG must be False.
        'SAVE_PLOTS': True,
        # Saves a plot of the evolution of the best fitness result for each
        # generation.

        # MULTIPROCESSING
        'MULTICORE': False,
        # Multiprocessing of phenotype evaluations.
        'CORES': cpu_count(),

        # STATE SAVING/LOADING
        'SAVE_STATE': False,
        # Saves the state of the evolutionary run every generation. You can
        # specify how often you want to save the state with SAVE_STATE_STEP.
        'SAVE_STATE_STEP': 1,
        # Specifies how often the state of the current evolutionary run is
        # saved (i.e. every n-th generation). Requires int value.
        'LOAD_STATE': None,
        # Loads an evolutionary run from a saved state. You must specify the
        # full file path to the desired state file. Note that state files have
        # no file type.

        # CACHING
        'CACHE': False,
        # The cache tracks unique individuals across evolution by saving a
        # string of each phenotype in a big list of all phenotypes. Saves all
        # fitness information on each individual. Gives you an idea of how much
        # repetition is in standard GE/GP.
        'LOOKUP_FITNESS': False,
        # Uses the cache to look up the fitness of duplicate individuals. CACHE
        #  must be set to True if you want to use this.
        'LOOKUP_BAD_FITNESS': False,
        # Uses the cache to give a bad fitness to duplicate individuals. CACHE
        # must be True if you want to use this (obviously)"""
        'MUTATE_DUPLICATES': False,
        # Removes duplicate individuals from the population by replacing them
        # with mutated versions of the original individual. Hopefully this will
        # encourage diversity in the population.

        # Set machine name (useful for doing multiple runs)
        'MACHINE': machine_name,

        # Set Random Seed
        'RANDOM_SEED': None
}


def load_params(file_name):
    """
    Load in a params text file and set the params dictionary directly.

    :param file_name: The name/location of a parameters file.
    :return: Nothing.
    """

    try:
        open(file_name, "r")
    except FileNotFoundError:
        print("Error: Parameters file not found. Ensure file\n"
              "       extension is specified, e.g. "
              "'regression.txt'.")
        quit()

    with open(file_name, 'r') as parameters:
        # Read the whole parameters file.
        content = parameters.readlines()

        for line in content:
            components = line.split(":")
            key, value = components[0], components[1].strip()

            # Evaluate parameters.
            try:
                value = eval(value)
            except:
                # We can't evaluate, leave value as a string.
                pass

            # Set parameter
            params[key] = value


def set_params(command_line_args):
    """
    This function parses all command line arguments specified by the user.
    If certain parameters are not set then defaults are used (e.g. random
    seeds, elite size). Sets the correct imports given command line
    arguments. Sets correct grammar file and fitness function. Also
    initialises save folders and tracker lists in utilities.trackers.

    :param command_line_args: Command line arguments specified by the user.
    :return: Nothing.
    """

    from utilities.algorithm.initialise_run import initialise_run_params
    from utilities.algorithm.initialise_run import set_param_imports
    from utilities.fitness.math_functions import return_percent
    from representation import grammar
    import utilities.algorithm.command_line_parser as parser
    from utilities.stats import trackers, clean_stats

    cmd_args, unknown = parser.parse_cmd_args(command_line_args)
    # TODO: how should we handle unknown parameters? Should we just add them indiscriminately to the params dictionary?

    # LOAD PARAMETERS FILE
    # NOTE that the parameters file overwrites all previously set parameters.
    if 'PARAMETERS' in cmd_args:
        load_params(path.join("..", "parameters", cmd_args['PARAMETERS']))

    # Join original params dictionary with command line specified arguments.
    # NOTE that command line arguments overwrite all previously set parameters.
    params.update(cmd_args)

    if params['LOAD_STATE']:
        # Load run from state.
        from utilities.algorithm.state import load_state
        
        # Load in state information.
        individuals = load_state(params['LOAD_STATE'])

        # Set correct search loop.
        from algorithm.search_loop import search_loop_from_state
        params['SEARCH_LOOP'] = search_loop_from_state
        
        # Set population.
        setattr(trackers, "state_individuals", individuals)
        
    else:
        # Elite size is set to either 1 or 1% of the population size, whichever is
        # bigger if no elite size is previously set.
        if params['ELITE_SIZE'] is None:
            params['ELITE_SIZE'] = return_percent(1, params['POPULATION_SIZE'])
    
        # Set the size of a generation
        params['GENERATION_SIZE'] = params['POPULATION_SIZE'] - params[
            'ELITE_SIZE']
    
        # Set GENOME_OPERATIONS automatically for faster linear operations.
        # TODO: there must be a cleaner way of doing this.
        if params['MUTATION'] in ['operators.mutation.int_flip', 'int_flip'] \
                and params['CROSSOVER'] in [
                    'operators.crossover.fixed_onepoint',
                    'operators.crossover.variable_onepoint',
                    'operators.crossover.fixed_twopoint',
                    'operators.crossover.variable_twopoint',
                    'fixed_onepoint', 'variable_onepoint',
                    'fixed_twopoint', 'variable_twopoint']:
            params['GENOME_OPERATIONS'] = True
        else:
            params['GENOME_OPERATIONS'] = False
            
        # Set correct param imports for specified function options, including
        # error metrics and fitness functions.
        set_param_imports()
    
        # Clean the stats dict to remove unused stats.
        clean_stats.clean_stats()
    
        # Initialise run lists and folders
        initialise_run_params()
    
        # Parse grammar file and set grammar class.
        params['BNF_GRAMMAR'] = grammar.Grammar(path.join("..", "grammars",
                                                params['GRAMMAR_FILE']))
