import matplotlib
import pandas as pd
from os import path, pathsep

from utilities.stats.trackers import best_fitness_list

matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rc('font', family='Times New Roman')


def save_best_fitness_plot():
    """
    Saves a plot of the current fitness.
    
    :return: Nothing
    """
    from algorithm.parameters import params

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(best_fitness_list)
    ax1.set_ylabel('fitness', fontsize=14)
    ax1.set_xlabel('Generation', fontsize=14)
    plt.title("Best fitness")
    plt.savefig(path.join(params['FILE_PATH'], "fitness.pdf"))
    plt.close()


def save_plot_from_data(data, name):
    """
    Saves a plot of a given set of data.
    
    :param data: the data to be plotted
    :param name: the name of the data to be plotted.
    :return: Nothing.
    """
    from algorithm.parameters import params

    # Plot the data
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(data)
    plt.title(name)
    plt.savefig(path.join(params['FILE_PATH'], (name + '.pdf')))
    plt.close()


def save_plot_from_file(filename, stat_name):
    """
    Saves a plot of a given stat from the stats file.
    
    :param filename: a full specified path to a .csv stats file.
    :param stat_name: the stat of interest for plotting.
    :return: Nothing.
    """

    # Read in the data
    data = pd.read_csv(filename, sep="\t")
    try:
        stat = list(data[stat_name])
    except KeyError:
        print("\nError: stat", stat_name, "does not exist")
        quit()

    # Plot the data
    fig = plt.figure()
    ax1 = fig.add_subplot(1, 1, 1)
    ax1.plot(stat)
    plt.title(stat_name)

    # Get save path
    save_path = pathsep.join(filename.split(pathsep)[:-1])

    # Save plot
    plt.savefig(path.join(save_path, (stat_name + '.pdf')))
    plt.close()
