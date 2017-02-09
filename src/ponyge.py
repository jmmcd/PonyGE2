#! /usr/bin/env python

# PonyGE
# Copyright (c) 2009 Erik Hemberg, James McDermott,
#                   Michael Fenton and David Fagan
# Hereby licensed under the GNU GPL v3.
""" Python GE implementation """

from utilities.algorithm.initialise_run import check_python_version

check_python_version()

from algorithm.parameters import params, set_params
from stats.stats import get_stats, stats
import sys


def mane():
    """ Run program """

    # Run evolution
    individuals = params['SEARCH_LOOP']()

    # Print final review
    get_stats(individuals, end=True)

    # Returns only needed if running experiment manager
    return params['TIME_STAMP'], stats['best_ever'].fitness


if __name__ == "__main__":
    set_params(sys.argv[1:]) # exclude the ponyge.py arg itself
    mane()
