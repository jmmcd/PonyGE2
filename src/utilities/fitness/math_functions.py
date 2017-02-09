import numpy as np
np.seterr(all="raise")


def return_percent(num, pop_size):
    """Returns either one percent of the population size or a given number,
       whichever is larger."""
    percent = int(round(pop_size/100))
    if percent < num:
        return num
    else:
        return percent


def aq(a, b):
    """The analytic quotient, intended as a "better protected division",
    from: Ji Ni and Russ H. Drieberg and Peter I. Rockett, "The Use of
    an Analytic Quotient Operator in Genetic Programming", IEEE
    Transactions on Evolutionary Computation."""
    return a / np.sqrt(1.0 + b**2.0)


def pdiv(x, y):
    # Koza's protected division
    # return 1 if y == 0 else x / y
    #
    # but we want this to work as an eval-able expression so we can't
    # include if-statements. Have to use some intermediate-level Numpy
    # trickery. This will raise a warning (depending on Numpy error
    # settings) because it always evaluates x / y before running
    # np.where.
    return np.where(y == 0, np.ones_like(x), x / y)


def rlog(x):
    # Koza's protected log: something not quite right...
    raise NotImplementedError
    # return np.where(x == 0, np.ones_like(x), np.log(np.abs(x)))


def ppow(x, y):
    # for the case where x negative and y non-integer
    return np.abs(x)**y


def ppow2(x, y):
    # for the case where x negative and y non-integer
    # preserve sign
    return np.sign(x) * (np.abs(x) ** y)


def psqrt(x):
    """ Protected square root operator"""
    return np.sqrt(np.abs(x))


def psqrt2(x):
    """ Protected square root operator"""
    # preserve sign
    return np.sign(x) * (np.sqrt(np.abs(x)))


def plog(x):
    """ Protected log operator"""
    return np.log(1.0 + np.abs(x))


def ave(x):
    """
    :param x: a given list
    :return: the average of param x
    """
    return np.mean(x)
