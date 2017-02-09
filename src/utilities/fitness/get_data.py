from os import listdir, getcwd, path

import numpy as np


def get_Xy_train_test(filename, randomise=True, test_proportion=0.5,
                      skip_header=0):
    """Read in a table of numbers and split it into X (all columns up
    to last) and y (last column), then split it into training and
    testing subsets according to test_proportion. Shuffle if
    required."""
    
    Xy = np.genfromtxt(filename, skip_header=skip_header)
    if randomise:
        np.random.shuffle(Xy)
    X = Xy[:, :-1]  # all columns but last
    y = Xy[:, -1]  # last column
    idx = int((1.0 - test_proportion) * len(y))
    train_X = X[:idx]
    train_y = y[:idx]
    test_X = X[idx:]
    test_y = y[idx:]
    
    return train_X, train_y, test_X, test_y


def get_Xy_train_test_separate(train_filename, test_filename, skip_header=0):
    """Read in training and testing data files, and split each into X
    (all columns up to last) and y (last column)."""

    # first try to auto-detect the field separator ie delimiter
    f = open(train_filename)
    for line in f:
        if line.startswith("#") or len(line) < 2:
            continue
        else:
            if "," in line:
                delimiter = ","; break
            elif "\t" in line:
                delimiter = "\t"; break
            elif ";" in line:
                delimiter = ";"; break
            elif ":" in line:
                delimiter = ":"; break
            else:
                delimiter = " "; break
    f.close()
    
    train_Xy = np.genfromtxt(train_filename, skip_header=skip_header,
                             delimiter=delimiter)
    test_Xy = np.genfromtxt(test_filename, skip_header=skip_header,
                            delimiter=delimiter)
    
    train_X = train_Xy[:, :-1].transpose()  # all columns but last
    train_y = train_Xy[:, -1].transpose()  # last column
    test_X = test_Xy[:, :-1].transpose()  # all columns but last
    test_y = test_Xy[:, -1].transpose()  # last column

    return train_X, train_y, test_X, test_y


def get_data(experiment, file_type="txt"):
    """ Return the training and test data for the current experiment.
    """

    datasets = listdir(getcwd() + "/../datasets/")
    for dataset in datasets:
        exp = dataset.split('.')[0].split('-')[0]
        if exp == experiment:
            file_type = dataset.split('.')[1]
    train_set = path.join("..", "datasets",
                          (experiment + "-Train." + str(file_type)))
    test_set = path.join("..", "datasets",
                         (experiment + "-Test." + str(file_type)))
    training_in, training_out, test_in, \
    test_out = get_Xy_train_test_separate(train_set, test_set, skip_header=1)
    
    return training_in, training_out, test_in, test_out
