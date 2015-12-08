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
from fnmatch import fnmatch
from collections import defaultdict
import csv
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("input", type=str,
                    help="A vcf file (version 4.0 or 4.1) or a folder of vcf files")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================


def snps_per_chromosome(infile):
    '''foo
    '''
    valid_chromosomes = set(["W","Z"] + [str(i) for i in range(1, 33)])
    CHROM = dict.fromkeys(valid_chromosomes, 0)
    with open(infile, "r") as vcf_reader:
        for line in vcf_reader:
            if line.startswith("#"):
                continue
            else:
                c = line.split()[0]
                if c in valid_chromosomes:
                    CHROM[c] += 1

    return CHROM


def handle_input(input):
    """If a directory is passed return all the files in the directory and sub-
     directories in a list, if a file is provided return a list with one entry.
    """
    if os.path.isdir(input):
        return list_files(input)
    elif os.path.isfile(input):
        return [input]


def split_tissue_infiles(infiles, tissues):

    tissue_files = {}
    for t in tissues:
        files = []
        for f in infiles:
            fname = os.path.split(f)[-1]
            if fname.startswith(t):
                files.append(f)
        tissue_files[t] = files

    return tissue_files


def combine_tissue_stats(infiles):

    filter_dict = {"coverage_filter": "",
                   "cluster_filter": "",
                   "known_filter": "",
                   "coding_filter": ""}

    coverage_pattern = "*_filteredallmin[0-9][0-9].vcf"
    cluster_pattern = "*_noclusters.vcf"
    known_pattern = "*_known.vcf"
    coding_pattern = "*_coding_transcript.vcf"

    # Depends on input file format
    for f in infiles:
        fname = os.path.split(f)[-1]
        if fnmatch(fname, coverage_pattern):
            filter_dict["coverage_filter"] = snps_per_chromosome(f)
        elif fnmatch(fname, cluster_pattern):
            filter_dict["cluster_filter"] = snps_per_chromosome(f)
        elif fnmatch(fname, known_pattern):
            filter_dict["known_filter"] = snps_per_chromosome(f)
        elif fnmatch(fname, coding_pattern):
            filter_dict["coding_filter"] = snps_per_chromosome(f)
        else:
            print "Could not match file: ", f

    return filter_dict


def write_chrom_stats(CHROM_STATS):
    valid_chromosomes = set(["W","Z"] + [str(i) for i in range(1, 33)])
    for tissue in CHROM_STATS:
        outfile = tissue+"-chromosome-snp-stats.csv"
        with open(outfile, "w") as f:
            chr_writer = csv.writer(f)
            header = ["coverage_filter", "cluster_filter", "known_filter", "coding_filter"]
            chr_writer.writerow(["CHR"]+header)
            for c in valid_chromosomes:
                row = [c]
                for h in header:
                    row.append(CHROM_STATS[tissue][h][c])
                chr_writer.writerow(row)



def main():

    # Prepare the function arguments for each infile
    fpattern = '*.vcf'
    infiles = [f for f in handle_input(args.input) if fnmatch(f, fpattern)]
    tissues = ["SPLEEN", "HEART", "GONAD", "LIVER"]
    tissue_files = split_tissue_infiles(infiles, tissues)

    CHROM_STATS = {}

    for t in tissue_files:
        CHROM_STATS[t] = combine_tissue_stats(tissue_files[t])

    write_chrom_stats(CHROM_STATS)


if __name__ == "__main__":
    main()
