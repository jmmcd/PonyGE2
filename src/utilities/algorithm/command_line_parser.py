import argparse


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
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""Welcome to PonyGE - Help
        -------------------
        The following are the available command line args
        please see src /algorithm/parameters.py
        for a more detailed explanation of each argument and possible
        values:""",
        epilog="""----------------------------
        To try out ponyge simply run: python ponyge.py

        Thanks for trying our product

        PonyGE Team""")

    # Set up class for checking float arguments
    class FloatAction(argparse.Action):
        """
        Class for checking a given float is within the range [0:1].
        """

        def __init__(self, option_strings, **kwargs):
            super(FloatAction, self).__init__(option_strings, **kwargs)

        def __call__(self, parser, namespace, value, option_string=None):
            if not 0 <= float(value) <= 1:
                print("\nError: parameter '", option_string,
                      "' outside allowed range [0:1]. "
                      "Value given:", value)
                quit()
            else:
                setattr(namespace, self.dest, float(value))

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

    # INDIVIDUAL SIZE
    parser.add_argument('--max_tree_depth',
                        dest='MAX_TREE_DEPTH',
                        type=int,
                        help='Sets the max derivation tree depth for the '
                             'algorithm, requires int value.')
    parser.add_argument('--max_tree_nodes',
                        dest='MAX_TREE_NODES',
                        type=int,
                        help='Sets the max derivation tree nodes for the '
                             'algorithm, requires int value.')
    parser.add_argument('--codon_size',
                        dest='CODON_SIZE',
                        type=int,
                        help='Sets the range from 0 to condon_size to be used '
                             'in genome, requires int value')
    parser.add_argument('--max_genome_length',
                        dest='MAX_GENOME_LENGTH',
                        type=int,
                        help='Sets the maximum chromosome length for the '
                             'algorithm, requires int value.')
    parser.add_argument('--max_wraps',
                        dest='MAX_WRAPS',
                        type=int,
                        help='Sets the maximum number of times the genome '
                             'mapping process can wrap over the length of the '
                             'genome. Requires int value.')

    # INITIALISATION
    parser.add_argument('--max_init_tree_depth',
                        dest='MAX_INIT_TREE_DEPTH',
                        type=int,
                        help='Sets the max tree depth for initialisation.')
    parser.add_argument('--min_init_tree_depth',
                        dest='MIN_INIT_TREE_DEPTH',
                        type=int,
                        help='Sets the min tree depth for initialisation.')
    parser.add_argument('--max_init_genome_length',
                        dest='MAX_INIT_GENOME_LENGTH',
                        type=int,
                        help='Sets the maximum length for chromosomes to be '
                             'initialised to, requires int value.')
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

    # EVALUATION
    parser.add_argument('--multicore',
                        dest='MULTICORE',
                        default=None,
                        action='store_true',
                        help='Turns on multicore evaluation.')
    parser.add_argument('--cores',
                        dest='CORES',
                        type=int,
                        help='Specify the number of cores to be used for '
                             'multicore evaluation. Requires int.')

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
                             'string such as "int_flip" or direct path string '
                             'such as "operators.mutation.int_flip".')
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
    parser.add_argument('--fitness_function',
                        dest='FITNESS_FUNCTION',
                        type=str,
                        help='Sets the fitness function to be used. '
                             'Requires string such as "regression".')
    parser.add_argument('--dataset',
                        dest='DATASET',
                        type=str,
                        help='For use with problems that use a dataset. '
                             'Requires string such as "Dow".')
    parser.add_argument('--target',
                        dest='TARGET',
                        type=str,
                        help='For string match problem. Requires target '
                             'string.')
    parser.add_argument('--experiment_name',
                        dest='EXPERIMENT_NAME',
                        type=str,
                        help='Optional parameter to save results in '
                             'results/[EXPERIMENT_NAME] folder. If not '
                             'specified then results are saved in default '
                             'results folder.')
    parser.add_argument('--error_metric',
                        dest='ERROR_METRIC',
                        type=str,
                        help='Sets the error metric to be used with regression'
                             ' style problems. Requires string such as "mse" '
                             'or "rmse".')
    parser.add_argument('--extra_fitness_parameters',
                        dest='EXTRA_FITNESS_PARAMETERS',
                        type=str,
                        help='Optional extra command line parameter for '
                             'inclusion of any extra information required '
                             'for user-specific runs. Can be whatever you '
                             'want it to be.')

    # OPTIONS
    parser.add_argument('--random_seed',
                        dest='RANDOM_SEED',
                        type=int,
                        help='Sets the seed to be used, requires int value.')
    parser.add_argument('--debug',
                        dest='DEBUG',
                        default=None,
                        action='store_true',
                        help='Disables saving of all ancillary files.')
    parser.add_argument('--verbose',
                        dest='VERBOSE',
                        default=None,
                        action='store_true',
                        help='Turns on the verbose output of the program in '
                             'terms of command line and extra files.')
    parser.add_argument('--silent',
                        dest='SILENT',
                        default=None,
                        action='store_true',
                        help='Prevents any output from being printed to the '
                             'command line.')
    parser.add_argument('--save_all',
                        dest='SAVE_ALL',
                        default=None,
                        action='store_true',
                        help='Saves the best phenotypes at each generation.')
    parser.add_argument('--save_plots',
                        dest='SAVE_PLOTS',
                        default=None,
                        action='store_true',
                        help='Saves plots for best fitness.')

    # STATE SAVING/LOADING
    parser.add_argument('--save_state',
                        dest='SAVE_STATE',
                        default=None,
                        action='store_true',
                        help='Saves the state of the evolutionary run every '
                             'generation. You can specify how often you want '
                             'to save the state with the command '
                             '"--save_state_step".')
    parser.add_argument('--save_state_step',
                        dest='SAVE_STATE_STEP',
                        default=None,
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
            setattr(namespace, 'LOOKUP_FITNESS', self.LOOKUP_FITNESS)
            setattr(namespace, 'LOOKUP_BAD_FITNESS', self.LOOKUP_BAD_FITNESS)
            setattr(namespace, 'MUTATE_DUPLICATES', self.MUTATE_DUPLICATES)

    # Generate a mutually exclusive group for caching options. This means
    # that you cannot specify multiple caching options simultaneously,
    # only one at a time.
    caching_group = parser.add_mutually_exclusive_group()
    caching_group.add_argument("--cache",
                               dest='CACHE',
                               action=CachingAction,
                               CACHE=True,
                               LOOKUP_FITNESS=True,
                               help='Tracks unique phenotypes and is used to '
                                    'lookup duplicate fitnesses.')
    caching_group.add_argument("--dont_lookup_fitness",
                               dest='CACHE',
                               action=CachingAction,
                               LOOKUP_FITNESS=False,
                               help='Turns on the cache to track duplicate '
                                    'individuals, but does not use the cache '
                                    'to save fitness evaluations.')
    caching_group.add_argument("--lookup_bad_fitness",
                               dest='CACHE',
                               action=CachingAction,
                               LOOKUP_FITNESS=True,
                               LOOKUP_BAD_FITNESS=True,
                               help='Gives duplicate phenotypes a bad fitness '
                                    'when encountered, requires cache.')
    caching_group.add_argument("--mutate_duplicates",
                               dest='CACHE',
                               action=CachingAction,
                               LOOKUP_FITNESS=False,
                               MUTATE_DUPLICATES=True,
                               help='Replaces duplicate individuals with '
                                    'mutated versions. Requires cache.')

    # Parse command line arguments using all above information.
    args, unknown = parser.parse_known_args(arguments)

    # All default args in the parser are set to "None". Only take arguments
    # which are not "None", i.e. arguments which have been passed in from
    # the command line.
    cmd_args = {key: value for key, value in vars(args).items() if value is
                not None}

    return cmd_args, unknown
