#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       run-filter-coverage.py
#==============================================================================
import argparse
import sys
import vcf
import os
# VCF requires Rpy2 internally and needs the correct R_HOME
os.environ["R_HOME"] = "/usr/lib64/R"
# Make utilities folder available
sys.path.append(os.path.abspath("../"))
from utilities.io import list_files
from utilities.runner import multiprocess
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("input", type=str,
                    help="A vcf file (version 4.0 or 4.1) or a folder of vcf files")
parser.add_argument("-c", "--coverage", type=int,
                    help="Minimum coverage of reads for a given position")
parser.add_argument("cov_table", type=str,
                    help="Variable coverage table based on binomial distribution")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================


def read_variable_coverage(cov_table):
    '''Read the variable coverage cutoff file computed using
       make-variable-coverage-table.py
    '''
    cov_cutoff = {}
    try:
        with open(cov_table, "rb") as infile:
            next(infile)  # Exclude header
            for line in infile:
                line = line.rstrip().split()
                cov_cutoff[int(line[0])] = int(line[1])
            return cov_cutoff
    except IOError:
        print "File does not exist: %s" % cov_table


def run_filter(arg_list):

    @multiprocess(arg_list)
    def filter_coverage(infile, outfile, cov_table, min_coverage=15):
        '''Filters VCF records by a constant minimum threshold (default 15) and
           by a variable threshold depending on the read depth. Currently a
           record is kept if any sample passes the filtering. For example, if
           three out of four samples do not pass the filtering, but one sample
           does the record will be written.
           This should guarantee, that at least one sample has the required
           coverage, while removing calls where every sample is invalid.

           2014-07-17: Added sanity check to ensure that all samples are
           present in the VCF file.
        '''
        vcf_reader = vcf.Reader(open(infile, "r"))
        vcf_writer = vcf.Writer(open(outfile, "w"), vcf_reader)

        max_cov = 0
        num_samples = len(vcf_reader.samples)
        missing_sample_records = 0

        for record in vcf_reader:

            # Ensure that all samples are included in the VCF file
            sample_count = (record.INFO["HET"] + record.INFO["NC"] +
                            record.INFO["HOM"] + record.INFO["WT"])

            if sample_count == num_samples:

                for sample in record.samples:

                    # Has the sample been called?
                    if sample["GT"] is not None:
                        num_ref = sample["RD"]
                        num_alt = sample["AD"]

                        coverage = num_ref + num_alt

                        if coverage not in cov_table:
                            continue

                        cutoff = cov_table[coverage]

                        if coverage >= min_coverage:
                            if num_ref >= cutoff and num_alt >= cutoff:
                                vcf_writer.write_record(record)
                                if coverage >= max_cov:
                                    max_cov = coverage
                                break
            else:
                missing_sample_records += 1

        print "Maximum coverage: ", max_cov
        print "Filtered records because of missing samples: ", missing_sample_records

    filter_coverage()


def handle_input(input):
    """If a directory is passed return all the files in the directory and sub-
     directories in a list, if a file is provided return a list with one entry.
    """
    if os.path.isdir(input):
        return list_files(input)
    elif os.path.isfile(input):
        return [input]


def main():
    cov_table = read_variable_coverage(args.cov_table)

    # Prepare the function arguments for each infile
    arg_list = []
    for infile in handle_input(args.input):

        head, tail = os.path.split(infile)
        outfile = tail.split(".")[0]+"_filteredmin"+str(args.coverage)+".vcf"
        outfile = os.path.join(head, outfile)
        print "Filtering file: %s" % infile

        arg_list.append((infile, outfile, cov_table, args.coverage))

    # Run the filtering step in parallel
    run_filter(arg_list)

if __name__ == "__main__":
    main()
