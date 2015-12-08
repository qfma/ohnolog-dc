#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       make-variable-coverage-table.py
#==============================================================================
import argparse
import sys
import scipy.stats
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("-c", "--coverage", type=int,
                    help="Maximum coverage of reads for a given position")
parser.add_argument("-e", "--error_rate", type=float,
                    help="Error rate for a read")
parser.add_argument("-t", "--threshold", type=float,
                    help="Probability threshold")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================


def calculate_variable_cutoff(max_cutoff, error_rate, threshold):
    '''
    The binomial distribution is actually pretty straightforward.
    If you have a coin that has a probability p to come up heads then
    he binomial distribution can tell us the probability that if you 
    throw n coins then what is the probability that you will get k heads.  
    We call this the probability mass function.
    Source: http://blog.mcbryan.co.uk/2013/02/the-binomial-distribution-python-and.html
    '''
    cutoffs = {}
    min_reads = 0
    for i in xrange(1, max_cutoff+1):
        for x in xrange(min_reads, i):
            prob = scipy.stats.binom.pmf(x, i, error_rate)
            if prob < threshold:
                min_reads = x
                cutoffs[i] = (x, prob)
                break
    return cutoffs


def calculate_variable_cutoff2(max_cutoff, error_rate, threshold):
    cutoffs = {}
    min_reads = 0
    for i in xrange(1, max_cutoff+1):
        for x in xrange(min_reads, i):
            pval = scipy.stats.binom_test(x, i, error_rate)
            if pval < threshold:
                min_reads = x
                cutoffs[i] = (x, pval)
                break
    return cutoffs


def main():
    # Minimum Phred score of 15 = 0.031622776601683791
    # 10**(-15./10.)
    cutoffs = calculate_variable_cutoff(args.coverage, args.error_rate, args.threshold)
    print "coverage", "minimum_reads", "pval"
    for c in cutoffs:
        print c, cutoffs[c][0], "{0:.15f}".format(cutoffs[c][1])

if __name__ == "__main__":
    main()
