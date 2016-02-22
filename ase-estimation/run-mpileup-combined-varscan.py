#!/usr/local/bin/python2.7
# -*- coding: utf-8 -*-
#
#       run-mpileup-combined-varscan.py
#==============================================================================
import argparse
import sys
import os
from subprocess import Popen, PIPE
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


def get_filtered_alignments(infolder):
    alns = []
    for root, dirs, files in os.walk(infolder):
        if "filter.bam" in files:
            alns.append(os.path.join(root, "filter.bam"))
    return alns


def get_snp_call_params(alns):

    genome = "../../FZDC-2014/rnaseq-mapping/data/ggal75/Gallus_gallus.Galgal4.75.dna.toplevel.fa"

    vcf = "1"
    min_var_feq = "1e-10"
    min_coverage = "2"
    min_qual = "20"
    pval = "1"
    strand = "0"
    depth = "10000000"
    # Set this to 90 for Z chromosome
    hom_freq = "0.90"
    # hom_freq = "0.75"

    samtools_params = ["samtools", "mpileup",
                       "-B",
                      # Does NOT need to be used, because mpileup outputs phred33
                      # "-6", # Assume the quality is in the Illumina-1.3+ encoding
                       "-d", depth,
                       # "-r", chromosome,
                       "-f", genome]

    # Add all the alignments to the mpileup command
    for a in alns:
        samtools_params.append(a)

    varscan_params = ["java", "-jar", "./tools/VarScan.v2.3.6.jar",
                      "mpileup2snp",
                      "--min-coverage", min_coverage,
                      "--min-avg-qual", min_qual,
                      "--strand-filter", strand,
                      "--p-value", pval,
                      "--min-freq-for-hom", hom_freq,
                      "--min-var-freq", min_var_feq,
                      "--output-vcf", vcf]
    return samtools_params, varscan_params


def run_snp_calls(alns, name):

    today = datetime.date.today()
    samtools, varscan = get_snp_call_params(alns)
    vcfname = str(today)+"_"+name+"_5mismatches.vcf"

    print vcfname
    with open(vcfname, "w") as vcfout:
        p1 = Popen(samtools, stdout=PIPE)
        p2 = Popen(varscan, stdin=p1.stdout, stdout=vcfout)
        p1.stdout.close()
        output = p2.communicate()[0]


def main():
    spleen_alns = get_filtered_alignments(args.spleen)
    gonad_alns = get_filtered_alignments(args.gonad)
    liver_alns = get_filtered_alignments(args.liver)
    heart_alns = get_filtered_alignments(args.heart)

    run_snp_calls(spleen_alns, "SPLEEN_M")
    run_snp_calls(gonad_alns, "GONAD_M")
    run_snp_calls(liver_alns, "LIVER_M")
    run_snp_calls(heart_alns, "HEART_M")

if __name__ == "__main__":
    main()
