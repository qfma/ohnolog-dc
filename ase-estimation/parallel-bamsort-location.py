#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#       parallel-bamsort-location.py
#==============================================================================
import argparse
import sys
import os
# Make utilities folder available
sys.path.append(os.path.abspath("../"))
from utilities.runner import exec_commands
from utilitis.io import list_files
#==============================================================================
#Command line options==========================================================
#==============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("indir", type=str,
                    help="A directory containing input BAM files")
if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()
#==============================================================================


def get_bam_files(files):
    bam_files = []
    for f in files:
        file_name, file_extension = os.path.splitext(f)
        if file_extension == ".bam":
            bam_files.append(f)
    return bam_files


def main():

    files = list_files(args.indir)
    bam_files = get_bam_files(files)
    commands = []

    for bam in bam_files:
        outfile = bam.split(".")[0]+".loc.sorted"
        commands.append(["samtools", "sort", bam, outfile])

    exec_commands(commands)

if __name__ == "__main__":
    main()
