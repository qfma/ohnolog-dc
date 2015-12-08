#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       combine-counts.py
#==============================================================================
import argparse
import sys
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import numpy as np
import datetime
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


def get_countfiles(infolder):
    counts = []
    for root, dirs, files in os.walk(infolder):
        for f in files:
            if f.endswith(".counts"):
                counts.append(os.path.join(root, f))
    return counts


def combine_counts(count_files):
    count_dataframes = []
    for f in count_files:
        sample_id = os.path.basename(f).split(".")[0]
        df = pd.read_csv(f, sep="\t", header=None, index_col=0)
        df.columns = [sample_id]
        # Convert everything to float, better for R
        df = df[[sample_id]].astype(float)
        # df.index.name = "Genes"
        count_dataframes.append(df)
    combined = pd.concat(count_dataframes,axis=1)
    # Drop the comments
    combined = combined.drop(["__no_feature", "__ambiguous", "__too_low_aQual",
                  "__not_aligned", "__alignment_not_unique"])

    # Reverse so that males are in front.
    # ID looks like this: WTCHG_35309_179_SPLEEN_M_1_sequence
    # Sorted by M and F
    sorted_samples = sorted(combined.columns, key = lambda x: x.split("_")[4], reverse=True)
    print sorted_samples
    combined = combined.reindex_axis(sorted_samples, axis=1)
    return combined


def main():
    spleen_counts = get_countfiles(args.spleen)
    heart_counts = get_countfiles(args.heart)
    liver_counts = get_countfiles(args.liver)
    gonad_counts = get_countfiles(args.gonad)

    spleen_combined_counts = combine_counts(spleen_counts)
    heart_combined_counts = combine_counts(heart_counts)
    liver_combined_counts = combine_counts(liver_counts)
    gonad_combined_counts = combine_counts(gonad_counts)




    today = datetime.date.today()
    spleen_combined_counts.to_csv(str(today)+"-MF-SPLEEN-COUNTS-COMBINED.expr")
    heart_combined_counts.to_csv(str(today)+"-MF-HEART-COUNTS-COMBINED.expr")
    liver_combined_counts.to_csv(str(today)+"-MF-LIVER-COUNTS-COMBINED.expr")
    gonad_combined_counts.to_csv(str(today)+"-MF-GONAD-COUNTS-COMBINED.expr")


    

if __name__ == "__main__":
    main()
