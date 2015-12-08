#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       index-bam.py
#==============================================================================
import argparse
import sys
import os
# Make utilities folder available
sys.path.append(os.path.abspath("../"))
from utilities.runner import exec_in_row
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("spleen", type=str,
                    help="A directory containing spleen files")
parser.add_argument("heart", type=str,
                    help="A directory containing heart files")
parser.add_argument("liver", type=str,
                    help="A directory containing liver files")
parser.add_argument("gonad", type=str,
                    help="A directory containing gonad files")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================


def get_filtered_alignments(infolder):
    alns = []
    for root, dirs, files in os.walk(infolder):
        if "filter.sorted.bam" in files:
            alns.append(os.path.join(root, "filter.sorted.bam"))
    return alns


def get_samtools_params(aln):
    samtools_params = ["samtools",
                    "index",
                    aln]
    return samtools_params


def run_samtools_index(all_alns):
    samtools_commands = []
    for aln in all_alns:
        samtools_commands.append(get_samtools_params(aln))
    
    exec_in_row(samtools_commands)

def main():
    spleen_alns = get_filtered_alignments(args.spleen)
    gonad_alns = get_filtered_alignments(args.gonad)
    liver_alns = get_filtered_alignments(args.liver)
    heart_alns = get_filtered_alignments(args.heart)

    run_samtools_index(spleen_alns)
    run_samtools_index(gonad_alns)
    run_samtools_index(liver_alns)
    run_samtools_index(heart_alns)

if __name__ == "__main__":
    main()
