#! /usr/bin/env python

# PonyGE2
# Copyright (c) 2017 Michael Fenton, James McDermott,
#                    David Fagan, Stefan Forstenlechner,
#                    and Erik Hemberg
# Hereby licensed under the GNU GPL v3.
""" Python GE implementation """

from PonyGE2.stats import get_stats
from PonyGE2.parameters import params

def mane():
    """ Run program """

    # Run evolution
    individuals = params['SEARCH_LOOP']()

    best = max(individuals)
    return best.phenotype, best.fitness

    # Print final review
    # get_stats(individuals, end=True)
