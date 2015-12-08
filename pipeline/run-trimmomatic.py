#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#       parallel-fastqc.py
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
        if f.endswith(".txt"):
            seqs.append(f)
    seqs.sort()
    pair1 = seqs[::2]
    pair2 = seqs[1::2]
    pairs = zip(pair1, pair2)

    return pairs


def run_trimmomatic(fastq_files):

    pairs = get_fastq_pairs(fastq_files)
    trimmomatic_commands = []

    for pair in pairs:
        trimmomatic_commands.append(trimmomatic_call(pair[0], pair[1]))

    print "Starting to run trimmomatic commands..."
    exec_in_row(trimmomatic_commands)


def trimmomatic_call(pair0, pair1):

    pair0_paired_trimmed = "%s.paired.trimmed" % pair0
    pair0_unpaired_trimmed = "%s.unpaired.trimmed" % pair0
    pair1_paired_trimmed = "%s.paired.trimmed" % pair1
    pair1_unpaired_trimmed = "%s.unpaired.trimmed" % pair1

    log = "%s.trimmed.log" % pair0
    leading = "LEADING:3"
    trailing = "TRAILING:3"
    swindow = "SLIDINGWINDOW:4:15"
    minlen = "MINLEN:36"
    trimmomatic_bin = "./tools/trimmomatic-0.22.jar"
    threads = "20"

    trimmomatic_params = ["java", "-classpath", trimmomatic_bin,
                          "org.usadellab.trimmomatic.TrimmomaticPE",
                          "-threads", threads,
                          "-trimlog", log,
                          pair0, pair1,
                          pair0_paired_trimmed, pair0_unpaired_trimmed,
                          pair1_paired_trimmed, pair1_unpaired_trimmed,
                          leading,
                          trailing,
                          swindow,
                          minlen]

    return trimmomatic_params

#==============================================================================


def main():
    spleen_files = list_folder(args.spleen)
    #heart_files = list_folder(args.heart)
    #liver_files = list_folder(args.liver)
    #gonad_files = list_folder(args.gonad)

    run_trimmomatic(spleen_files)
    #run_trimmomatic(heart_files)
    #run_trimmomatic(liver_files)
    #run_trimmomatic(gonad_files)

if __name__ == "__main__":
    main()
