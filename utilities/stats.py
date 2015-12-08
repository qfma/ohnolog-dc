#!/usr/bin/env python2.7
''' Mathematical functions

This file provides short mathematical functions that have not been
implemented in the standard library.

'''

from math import pow
from math import log


def reverse_log2(l):
    '''Unlog every value in list'''
    return [pow(2, i) for i in l]


def log2(i):
    '''Log2 transform the expression, values below 2 RPKM are set
    to 0 to avoid negative expression values und -inf errors.
    This is maybe a really bad idea!'''
    if i == 0:
        return float(0.0)
    elif i > 0:
        return float(log(i, 2))
