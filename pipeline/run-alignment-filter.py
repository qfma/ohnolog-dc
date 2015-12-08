#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       run-tophat.py
#==============================================================================
import argparse
import sys
import pysam
import os
from utilities.runner import multiprocess
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


def get_alignments(infolder):
    alns = []
    for root, dirs, files in os.walk(infolder):
        if "accepted_hits.bam" in files:
            alns.append(os.path.join(root, "accepted_hits.bam"))
    return alns


def run_alignment_filter(alns):

    @multiprocess(alns)
    def remove_non_unique_reads(aln):

        head, tail = os.path.split(aln)
        filter_file = os.path.join(head, "filter.bam")

        raw_reads = pysam.Samfile(aln, 'rb')
        filter_reads = pysam.Samfile(filter_file, 'wb', template=raw_reads)

        print "Filtering %s..." % aln
        for read in raw_reads.fetch():
            if ('NH', 1) in read.tags:
                filter_reads.write(read)

        raw_reads.close()
        filter_reads.close()
        pysam.index(filter_file)

    @multiprocess(alns)
    def index_aln(aln):

        head, tail = os.path.split(aln)
        bai_file = os.path.join(head, "accepted_hits.bam.bai")

        if not os.path.exists(bai_file):
            print "Indexing %s..." % aln
            pysam.index(aln)
        else:
            print "BAM index exists already...continue"

    index_aln()
    remove_non_unique_reads()


def main():
    spleen_alns = get_alignments(args.spleen)
    gonad_alns = get_alignments(args.gonad)
    liver_alns = get_alignments(args.liver)
    heart_alns = get_alignments(args.heart)

    all_alns = spleen_alns + gonad_alns + liver_alns + heart_alns
    run_alignment_filter(all_alns)

if __name__ == "__main__":
    main()
