import getopt
import sys
from datetime import datetime
from os import getcwd, listdir, path, pathsep

import matplotlib
import numpy as np
np.seterr(all="raise")
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('font', family='Times New Roman')


def help_message():
    """
    Prints a help message to explain the usage of this file.
    
    :return: Nothing
    """
    
    lines_1 = ["Welcome to PonyGE's post-run stats parser.",
               "-------------------",
               "The following are the available command line args. You must "
               "specify an experiment name and at least one stat to be parsed."
               ""]
    
    lines_2 = [["\t--help:", "Shows this help message."],
               ["\t--experiment_name:", "The name of the containing folder in "
                                        "which target runs are saved, e.g. "
                                        "in /results/[EXPERIMENT_NAME]."],
               ["\t--stats:", "A comma-separated list of stats to be parsed. "
                              "Specified stats must be valid keys from the "
                              "stats.stats.stats dictionary. "
                              "IMPORTANT: MUST NOT CONTAIN ANY SPACES."],
               ["\t--graph:", "Saves a .pdf figure of each stat specified."]]

    # This simply justifies the print statement such that it is visually
    # pleasing to look at.
    for line in lines_1:
        print(line)
    col_width = max(len(line[0]) for line in lines_2)
    for line in sorted(lines_2):
        print(" ".join(words.ljust(col_width) for words in line))


def parse_opts(command_line_args):
    """
    Parses command line arguments and returns usable variables which are used
    as inputs for other functions in this file.
    
    :param command_line_args: flags passed in from the command line at
    execution
    :return: experiment_name: the name of the containing folder of results,
             stats: a list of strings. Each string is a valid key from the
                    stats.stats.stats dictionary.
             graph: an optional boolean flag for graphing specified stats.
    """
    
    try:
        opts, args = getopt.getopt(command_line_args[1:], "",
                                   ["help", "experiment_name=", "stats=",
                                    "graph"])
    except getopt.GetoptError as err:
        print("In order to parse stats you need to specify the location of the"
              " target stats files and a list of desired stats to parse. \n",
              "Run python parse_stats.py --help for more info")
        print(str(err))
        exit(2)

    if not opts:
        print("In order to parse stats you need to specify the location of the"
              " target stats files and a list of desired stats to parse. \n",
              "Run python parse_stats.py --help for more info")
        exit(2)

    experiment_name, stats, graph = None, None, False
    
    # iterate over all arguments in the option parser.
    for opt, arg in opts:
        if opt == "--help":
            # Print help message.
            help_message()
            exit()
        
        elif opt == "--experiment_name":
            # Set experiment name (i.e. containing folder for multiple runs).
            experiment_name = arg
        
        elif opt == "--stats":
            # List stats. Stats must be parsed correctly.
            # TODO: There must be a better way to pass in a list of options.
            if arg[0] == "[" and arg[-1] == "]":
                stats = arg[1:-1].split(",")
            
            else:
                stats = arg.split(",")
        
        elif opt == "--graph":
            # Set boolean flag for graphing stats.
            graph = True
    
    return experiment_name, stats, graph


def parse_stat_from_runs(experiment_name, stats, graph):
    """
    Analyses a list of given stats from a group of runs saved under an
    "experiment_name" folder. Creates a summary .csv file which can be used by
    plotting functions in utilities.save_plot. Saves a file of the format:

        run0_gen0       run1_gen0       .   .   .   run(n-1)_gen0
        run0_gen1       run1_gen1       .   .   .   run(n-1)_gen1
        run0_gen2       run1_gen2       .   .   .   run(n-1)_gen2
        .               .               .   .   .   .
        .               .               .   .   .   .
        .               .               .   .   .   .
        run0_gen(n-1)   run1_gen(n-1)   .   .   .   run(n-1)_gen(n-1)
        run0_gen(n)     run1_gen(n)     .   .   .   run(n-1)_gen(n)
        
    Generated file is compatible with
        
        utilities.save_plot.save_average_plot_across_runs()
    
    :param experiment_name: The name of a collecting folder within the
    ./results folder which holds multiple runs.
    :param stats: A list of the names of the stats to be parsed.
    :param graph: A boolean flag for whether or not to save figure.
    :return: Nothing.
    """

    # Since results files are not kept in source directory, need to escape
    # one folder.
    file_path = path.join(getcwd(), "..", "results")

    # Check for use of experiment manager.
    if experiment_name:
        file_path = path.join(file_path, experiment_name)
    
    else:
        print("Error: experiment name not specified")
        quit()
    
    # Find list of all runs contained in the specified folder.
    runs = [run for run in listdir(file_path) if "." not in run]

    # Place to store the header for full stats file.
    header = ""

    # Array to store all stats
    full_stats = []

    for stat in stats:
        # Iterate over all specified stats.

        if stat == "best_ever":
            print("Error: Cannot graph instances of individual class. Do not"
                  " specify 'best_ever' as stat to be parsed.")
            quit()

        summary_stats = []

        # Iterate over all runs
        for run in runs:
            # Get file name
            file_name = path.join(file_path, str(run), "stats.tsv")

            # Load in data
            data = pd.read_csv(file_name, sep="\t")

            try:
                # Try to extract specific stat from the data.
                summary_stats.append(list(data[stat]))

            except KeyError:
                # The requested stat doesn't exist.
                print("Error: stat", stat, "does not exist")
                quit()

        if stat in ["total_time", "time_taken"]:
            # Must parse time-related stats to the correct time format.
            zero = datetime.strptime("1900-01-01 0:00:00.000000",
                                     "%Y-%m-%d %H:%M:%S.%f")
            for i, run in enumerate(summary_stats):
                summary_stats[i] = [(datetime.strptime(time,
                                                       "%H:%M:%S.%f") -
                                     zero).total_seconds()
                                    for time in run]

        # Generate numpy array of all stats
        summary_stats = np.array(summary_stats)

        # Append Stat to header.
        header = header + stat + "_mean" + ","

        summary_stats_mean = np.nanmean(summary_stats, axis=0)
        full_stats.append(summary_stats_mean)

        # Append Stat to header.
        header = header + stat + "_std" + ","
        summary_stats_std = np.nanstd(summary_stats, axis=0)
        full_stats.append(summary_stats_std)
        summary_stats = np.transpose(summary_stats)

        # Save stats as a .csv file.
        np.savetxt(path.join(file_path, (stat + ".csv")), summary_stats,
                   delimiter=",")

        if graph:
            # Graph stat by calling graphing function.
            save_average_plot_across_runs(path.join(file_path, (stat +
                                                                ".csv")))
    # Convert and rotate full stats
    full_stats = np.array(full_stats)
    full_stats = np.transpose(full_stats)

    # Save full stats to csv file.
    np.savetxt(path.join(file_path, "full_stats.csv"), full_stats,
               delimiter = ",", header=header[:-1])


def save_average_plot_across_runs(filename):
    """
    Saves an average plot of multiple runs. Input file data must be of the
    format:

        run0_gen0       run1_gen0       .   .   .   run(n-1)_gen0
        run0_gen1       run1_gen1       .   .   .   run(n-1)_gen1
        run0_gen2       run1_gen2       .   .   .   run(n-1)_gen2
        .               .               .   .   .   .
        .               .               .   .   .   .
        .               .               .   .   .   .
        run0_gen(n-1)   run1_gen(n-1)   .   .   .   run(n-1)_gen(n-1)
        run0_gen(n)     run1_gen(n)     .   .   .   run(n-1)_gen(n)

    The required file can be generated using

        stats.parse_stats.parse_stat_from_runs()

    Generates a .pdf graph of average value with standard deviation.
    :param filename: the full file name of a .csv file containing the fitnesses
    of each generation of multiple runs. Must be comma separated.
    :return: Nothing.
    """
    
    # Get stat name from the filename. Used later for saving graph.
    stat_name = filename.split(pathsep)[-1].split(".")[0]
    
    # Load in data.
    data = np.genfromtxt(filename, delimiter=',')[:, :-1]
    
    # Generate average and standard deviations of loaded data.
    ave = np.nanmean(data, axis=1)
    std = np.nanstd(data, axis=1)
    
    # Calculate max and min of standard deviation.
    stdmax = ave + std
    stdmin = ave - std
    
    # Generate generation range over which data is to be graphed.
    max_gens = len(ave)
    r = range(1, max_gens + 1)
    
    # Initialise figure plot.
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    
    # Plot data and standard deviation infill.
    ax1.plot(r, ave, color="blue")
    ax1.fill_between(r, stdmin, stdmax, color="DodgerBlue", alpha=0.5)

    # Set x-axis limits.
    plt.xlim(0, max_gens + 1)
    
    # Set title and axes.
    plt.title("Average " + stat_name)
    plt.xlabel('Generation', fontsize=14)
    plt.ylabel('Average ' + stat_name, fontsize=14)
    
    # Save graph under the same name as the original .csv file but with a
    # .pdf extension instead.
    new_filename = filename[:-3] + "pdf"
    plt.savefig(str(new_filename))
    
    plt.close()


if __name__ == "__main__":
    experiment_name, stats, graph = parse_opts(sys.argv)
    parse_stat_from_runs(experiment_name, stats, graph)
