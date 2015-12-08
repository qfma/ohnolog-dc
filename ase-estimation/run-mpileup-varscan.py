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
        if "filter.bam" in files:
            alns.append(os.path.join(root, "filter.bam"))
    return alns


def get_snp_call_params(aln):

    genome = "../rnaseq-pipeline/ggal75/Gallus_gallus.Galgal4.75.dna.toplevel.fa"

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
                       "-d", depth,
                       # "-r", chromosome,
                       "-f", genome,
                       aln]

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


def run_snp_calls(all_alns):

    @multiprocess(all_alns)
    def snp_call(aln):
        sample_id = aln.split("/")[-2]
        samtools, varscan = get_snp_call_params(aln)
        vcfname = sample_id+"_5mismatches.vcf"
        print vcfname
        with open(vcfname, "wb") as vcfout:
            p1 = Popen(samtools, stdout=PIPE)
            p2 = Popen(varscan, stdin=p1.stdout, stdout=vcfout)
            p1.stdout.close()
            output = p2.communicate()[0]

    snp_call()


def main():
    spleen_alns = get_filtered_alignments(args.spleen)
    gonad_alns = get_filtered_alignments(args.gonad)
    liver_alns = get_filtered_alignments(args.liver)
    heart_alns = get_filtered_alignments(args.heart)

    run_snp_calls(spleen_alns)
    run_snp_calls(gonad_alns)
    run_snp_calls(liver_alns)
    run_snp_calls(heart_alns)

if __name__ == "__main__":
    main()
