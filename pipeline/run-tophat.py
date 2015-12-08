#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       run-tophat.py
#==============================================================================
import argparse
import sys
# sys.path.append(os.path.abspath(".."))
from utilities.runner import exec_in_row
from utilities.io import list_folder
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


def get_fastq_pairs(fastq_files):
    seqs = []
    for f in fastq_files:
        if f.split(".")[-1] == "trimmed" and f.split(".")[-2] == "paired":
            seqs.append(f)
    seqs.sort()

    pair1 = seqs[::2]
    pair2 = seqs[1::2]
    pairs = zip(pair1, pair2)

    return pairs


def run_spleen(fastq_files):

    pairs = get_fastq_pairs(fastq_files)
    tophat_commands = []

    for pair in pairs:
        mate0 = pair[0]
        mate1 = pair[1]
        basename = mate0.split(".")[0]
        output = basename
        # if not os.path.isdir(output):
        #     os.mkdir(output)
        tophat_commands.append(tophat_call(output, mate0, mate1))

    print "Starting to run tophat2 commands..."
    print tophat_commands
    exec_in_row(tophat_commands)


def run_heart(fastq_files):

    pairs = get_fastq_pairs(fastq_files)
    tophat_commands = []

    for pair in pairs:
        mate0 = pair[0]
        mate1 = pair[1]
        basename = mate0.split(".")[0]
        output = basename
        # if not os.path.isdir(output):
        #     os.mkdir(output)
        tophat_commands.append(tophat_call(output, mate0, mate1))

    print "Starting to run tophat2 commands..."
    print tophat_commands
    exec_in_row(tophat_commands)


def run_liver(fastq_files):

    pairs = get_fastq_pairs(fastq_files)
    tophat_commands = []

    for pair in pairs:
        mate0 = pair[0]
        mate1 = pair[1]
        basename = mate0.split(".")[0]
        output = basename
        # if not os.path.isdir(output):
        #     os.mkdir(output)
        tophat_commands.append(tophat_call(output, mate0, mate1))

    print "Starting to run tophat2 commands..."
    print tophat_commands
    exec_in_row(tophat_commands)


def run_gonad(fastq_files):

    pairs = get_fastq_pairs(fastq_files)
    tophat_commands = []

    for pair in pairs:
        mate0 = pair[0]
        mate1 = pair[1]
        basename = mate0.split(".")[0]
        output = basename
        # if not os.path.isdir(output):
        #     os.mkdir(output)
        tophat_commands.append(tophat_call(output, mate0, mate1))

    print "Starting to run tophat2 commands..."
    print tophat_commands
    exec_in_row(tophat_commands)


def tophat_call(output, mate0, mate1):

    gtf = "./data/ggal75/Gallus_gallus.Galgal4.75.NORNA.gtf"
    genome = "./data/ggal75/Gallus_gallus.Galgal4.75.dna.toplevel"
    mismatches = "5"
    threads = "22"

    tophat_params = ["./tools/tophat-2.0.11.Linux_x86_64/tophat2",
                     "-o", output,
                     "--solexa1.3-quals",
                     "-p", threads,
                     "-G", gtf,
                     "-N", mismatches,
                     "--b2-L", "20",
                     "--b2-N", "1",
                     "--read-edit-dist", mismatches,
                     "--no-novel-juncs",
                     genome,
                     mate0,
                     mate1]

    return tophat_params


def main():
    spleen_files = list_folder(args.spleen)
    heart_files = list_folder(args.heart)
    liver_files = list_folder(args.liver)
    gonad_files = list_folder(args.gonad)

 #   run_spleen(spleen_files)
    run_heart(heart_files)
    run_liver(liver_files)
    run_gonad(gonad_files)

if __name__ == "__main__":
    main()
