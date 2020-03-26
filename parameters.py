"""Algorithm parameters"""
from multiprocessing import cpu_count
from socket import gethostname

hostname = gethostname().split('.')
machine_name = hostname[0]

params = {
        # Set default step and search loop functions
        'SEARCH_LOOP': 'PonyGE2.algorithm.search_loop',
        'STEP': 'PonyGE2.algorithm.step',

        # Evolutionary Parameters
        'POPULATION_SIZE': 10,
        'GENERATIONS': 500,
        'HILL_CLIMBING_HISTORY': 1000,
        'SCHC_COUNT_METHOD': "count_all",

        # Set optional experiment name
        'EXPERIMENT_NAME': None,
        # Set default number of runs to be done.
        # ONLY USED WITH EXPERIMENT MANAGER.
        'RUNS': 1,

        # Class of problem
        'FITNESS_FUNCTION': "arc_evaluate",

        # Select problem dataset
        'DATASET_TRAIN': "../../training/25ff71a9.json",
        'DATASET_TEST': "../../test/25ff71a9.json",
        'DATASET_DELIMITER': None,

        # Set grammar file
        'GRAMMAR_FILE': "PonyGE2/grammars/arc.bnf",

        # Set the number of depths permutations are calculated for
        # (starting from the minimum path of the grammar).
        # Mainly for use with the grammar analyser script.
        'PERMUTATION_RAMPS': 5,

        # Select error metric
        'ERROR_METRIC': None,

        # Optimise constants in the supervised_learning fitness function.
        'OPTIMIZE_CONSTANTS': False,

        # Specify target for target problems
        'TARGET': "ponyge_rocks",

        # Set max sizes of individuals
        'MAX_TREE_DEPTH': 30,  # SET TO 90 DUE TO PYTHON EVAL() STACK LIMIT.
                               # INCREASE AT YOUR OWN RISK.
        'MAX_TREE_NODES': None,
        'CODON_SIZE': 1000,
        'MAX_GENOME_LENGTH': 500,
        'MAX_WRAPS': 0,

        # INITIALISATION
        # Set initialisation operator.
        'INITIALISATION': "PonyGE2.operators.uniform_tree",
        # Set the maximum geneome length for initialisation.
        'INIT_GENOME_LENGTH': 200,
        # Set the maximum tree depth for initialisation.
        'MAX_INIT_TREE_DEPTH': 5,
        # Set the minimum tree depth for initialisation.
        'MIN_INIT_TREE_DEPTH': None,

        # SELECTION
        # Set selection operator.
        'SELECTION': "PonyGE2.operators.tournament",
        # For tournament selection
        'TOURNAMENT_SIZE': 2,
        # For truncation selection
        'SELECTION_PROPORTION': 0.5,
        # Allow for selection of invalid individuals during selection process.
        'INVALID_SELECTION': False,

        # OPERATOR OPTIONS
        # Boolean flag for selecting whether or not mutation is confined to
        # within the used portion of the genome. Default set to True.
        'WITHIN_USED': True,

        # CROSSOVER
        # Set crossover operator.
        'CROSSOVER': "PonyGE2.operators.variable_onepoint",
        # Set crossover probability.
        'CROSSOVER_PROBABILITY': 0.5,
        # Prevents crossover from generating invalids.
        'NO_CROSSOVER_INVALIDS': False,

        # MUTATION
        # Set mutation operator.
        'MUTATION': "PonyGE2.operators.int_flip_per_codon",
        # Set mutation probability (None defaults to 1 over the length of
        # the genome for each codon)
        'MUTATION_PROBABILITY': 0.25,
        # Set number of mutation events
        'MUTATION_EVENTS': 1,
        # Prevents mutation from generating invalids.
        'NO_MUTATION_INVALIDS': False,

        # REPLACEMENT
        # Set replacement operator.
        'REPLACEMENT': "PonyGE2.operators.generational",
        # Set elite size.
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
        # Save the phenotype of the best individual from each generation. Can
        # generate a lot of files. DEBUG must be False.
        'SAVE_ALL': False,
        # Save a plot of the evolution of the best fitness result for each
        # generation.
        'SAVE_PLOTS': False,

        # MULTIPROCESSING
        # Multi-core parallel processing of phenotype evaluations.
        'MULTICORE': True,
        # Set the number of cpus to be used for multiprocessing
        'CORES': cpu_count(),

        # STATE SAVING/LOADING
        # Save the state of the evolutionary run every generation. You can
        # specify how often you want to save the state with SAVE_STATE_STEP.
        'SAVE_STATE': False,
        # Specify how often the state of the current evolutionary run is
        # saved (i.e. every n-th generation). Requires int value.
        'SAVE_STATE_STEP': 1,
        # Load an evolutionary run from a saved state. You must specify the
        # full file path to the desired state file. Note that state files have
        # no file type.
        'LOAD_STATE': None,

        # SEEDING
        # Specify a list of PonyGE2 individuals with which to seed the initial
        # population.
        'SEED_INDIVIDUALS': [],
        # Specify a target seed folder in the 'seeds' directory that contains a
        # population of individuals with which to seed a run.
        'TARGET_SEED_FOLDER': None,
        # Set a target phenotype string for reverse mapping into a GE
        # individual
        'REVERSE_MAPPING_TARGET': None,
        # Set Random Seed for all Random Number Generators to be used by
        # PonyGE2, including the standard Python RNG and the NumPy RNG.
        'RANDOM_SEED': None,

        # CACHING
        # The cache tracks unique individuals across evolution by saving a
        # string of each phenotype in a big list of all phenotypes. Saves all
        # fitness information on each individual. Gives you an idea of how much
        # repetition is in standard GE/GP.
        'CACHE': True,
        # Uses the cache to look up the fitness of duplicate individuals. CACHE
        # must be set to True if you want to use this.
        'LOOKUP_FITNESS': False,
        # Uses the cache to give a bad fitness to duplicate individuals. CACHE
        # must be True if you want to use this (obviously)
        'LOOKUP_BAD_FITNESS': False,
        # Removes duplicate individuals from the population by replacing them
        # with mutated versions of the original individual. Hopefully this will
        # encourage diversity in the population.
        'MUTATE_DUPLICATES': False,

        # MULTIAGENT Parameters
        # True or False for Multiagent
        'MULTIAGENT': False,
        # Agent Size. Number of agents having their own copy of genetic material
        'AGENT_SIZE': 100,
        # Interaction Probablity. How frequently the agents can interaction with each other
        'INTERACTION_PROBABILITY': 0.5,
        
        # OTHER
        # Set machine name (useful for doing multiple runs)
        'MACHINE': machine_name
}