#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       run-mpileup-varscan.py
#==============================================================================
import argparse
import sys
import os
# Make utilities folder available
sys.path.append(os.path.abspath("../"))
from utilities.runner import multiprocess
from subprocess import Popen, PIPE
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


def get_htseq_params(aln):

    gtf = "../../FZDC-2014/rnaseq-mapping/data/ggal75/Gallus_gallus.Galgal4.75.NORNA.gtf"

    htseq_params = ["htseq-count",
                    "-f", "bam",  # Use bam files
                    "-r", "name",  # Tophat sorts by position! Needs to be sorted by name using samtools!
                    "-s", "no",   # We do not have stranded data
                    aln, gtf]
    return htseq_params


def run_count_extraction(all_alns):

    @multiprocess(all_alns)
    def count_extraction(aln):

        sample_id = aln.split("/")[-2]
        print sample_id
        htseq = get_htseq_params(aln)
        print htseq
        print aln
        countname = sample_id+".counts"
        errname = sample_id+".htseqerr"
        print countname, errname
        print "---"
        with open(countname, "wb") as countout, open(errname, "wb") as counterr:
            p1 = Popen(htseq,
                       stdout=countout,  # Output file name
                       stderr=counterr   # Htseq error output
                       )

    count_extraction()


def main():
    spleen_alns = get_filtered_alignments(args.spleen)
    gonad_alns = get_filtered_alignments(args.gonad)
    liver_alns = get_filtered_alignments(args.liver)
    heart_alns = get_filtered_alignments(args.heart)

    run_count_extraction(spleen_alns)
    run_count_extraction(gonad_alns)
    run_count_extraction(liver_alns)
    run_count_extraction(heart_alns)

if __name__ == "__main__":
    main()
