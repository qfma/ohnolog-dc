#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       extract-feature-from-vcf.py
#==============================================================================
import argparse
import sys
import os
# Make utilities folder available
sys.path.append(os.path.abspath("../"))
from utilities.io import list_files
from utilities.runner import multiprocess
from fnmatch import fnmatch
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("input", type=str,
                    help="A vcf file (version 4.0 or 4.1) or a folder of vcf files")
parser.add_argument("-f", "--feature", type=str,
                    help="The feature to be extracted eg. gene, transcript")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================


def run_feature_extract(arg_list):

    @multiprocess(arg_list)
    def feature_extract(infile, outfile, feature):
        '''foo
        '''
        with open(infile, "r") as vcf_reader, open(outfile, "w") as vcf_writer:
            for line in vcf_reader:
                if line.startswith("##"):
                    vcf_writer.write(line)
                elif feature in line.split():
                    vcf_writer.write(line)

    feature_extract()


def handle_input(input):
    """If a directory is passed return all the files in the directory and sub-
     directories in a list, if a file is provided return a list with one entry.
    """
    if os.path.isdir(input):
        return list_files(input)
    elif os.path.isfile(input):
        return [input]


def main():

    # Prepare the function arguments for each infile
    arg_list = []
    fpattern = '*_coding.vcf'
    infiles = [f for f in handle_input(args.input) if fnmatch(f, fpattern)]
    for infile in infiles:

        head, tail = os.path.split(infile)
        outfile = tail.split(".")[0]+"_"+str(args.feature)+".vcf"
        outfile = os.path.join(head, outfile)
        print "Filtering file: %s" % infile

        arg_list.append((infile, outfile, args.feature))

    # Run the filtering step in parallel
    run_feature_extract(arg_list)

if __name__ == "__main__":
    main()
