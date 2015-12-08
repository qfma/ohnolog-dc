#!/usr/bin/env python2.7
''' Miscellaneous helpers

I cannot come up with a good group at the moment..
'''
from collections import Counter
import operator


def avg(l):
    return sum(l) / float(len(l))


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def count_elements(l):
    '''Count the number of reoccuring elements in a list'''
    return Counter(l)


def get_dict_max(d):
    '''Returns the key of the item in
       d with the highest value'''
    return max(d, key=d.get)


def get_dict_min(d):
    '''Returns the key of the item in
       d with the lowest value'''
    return min(d, key=d.get)


def chunks(l, n):
    """ Yield successive n-sized chunks from l.
    """
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def flatten(l):
    return [item for sublist in l for item in sublist]
