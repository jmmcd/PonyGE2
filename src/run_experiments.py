""" This program cues up and executes multiple runs of PYGE. Results of runs
    are parsed and placed in a spreadsheet for easy visual analysis.

    Copyright (c) 2014 Michael Fenton
    Hereby licensed under the GNU GPL v3."""

from utilities.algorithm.initialise_run import check_python_version

check_python_version()

from multiprocessing import cpu_count, Pool
from operator import itemgetter
from datetime import datetime
from subprocess import call
import sys

from algorithm.parameters import params, set_params
from stats.parse_stats import parse_stat_from_runs


def execute_run(seed):
    """
    Initialise all aspects of a run.
    
    :return: Nothing.
    """

    exec_str = "source activate py35;python ponyge.py " \
               "--random_seed " + str(seed) + " ".join(sys.argv[1:])
    call(exec_str, shell=True)


def execute_multi_single_core():
    """ Execute multiple runs in series using multiple cores to evaluate each
        population. """

    results, run_results = [], []
    time1 = datetime.now()
    pool = Pool(processes=CORES)
    
    print("\nMulti-Run Start:", time1, "\n")

    for run in range(RUNS):

        # Execute a single evolutionary run.
        results.append(pool.apply_async(execute_run, (run)))
    
    for result in results:
        result.get()
    
    pool.close()
    
    # Save spreadsheet and average fitness plot for runs so far.
    parse_stat_from_runs(params['EXPERIMENT_NAME'],
                         ["best_fitness"], True)

    time2 = datetime.now()
    total_time = time2 - time1

    # Write info about best indiv from each run to a file.
    filename = "../results/" + str(params['EXPERIMENT_NAME']) + "/" + \
               str(params['EXPERIMENT_NAME']) + ".txt"
    savefile = open(filename, 'w')
    for ans, answer in enumerate(run_results):
        savefile.write("Run " + str(ans) + "\tName: " + str(answer[0]) +
                       "\tBest: " + str(answer[1]) + "\n")
    run_results.sort(key=itemgetter(1))

    # Write info about best overall result.
    savefile.write("\nBEST: " + str(run_results[-1]))
    savefile.write("\n\nTotal time taken for " + str(RUNS) +
                   " runs: " + str(total_time))
    savefile.close()

    print("\nTotal time taken for", RUNS, "runs:", total_time)
    print("\nBEST:", run_results[-1])


if __name__ == "__main__":
    # Setup run parameters.
    set_params(sys.argv)
    RUNS = 30
    CORES = cpu_count() - 1

    # Execute multiple runs.
    execute_multi_single_core()
