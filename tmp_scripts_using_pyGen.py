#! /usr/bin/env python

# IMPORT
#-------------------------------------------------------------------------------
import argparse
import logging as log


# FUNCTIONS
#-------------------------------------------------------------------------------
def snpAnalyse(gtf, vcf, gos, snpAnnotOut):
    import SNPs_annotator as snp

    log.info("STARTING snps pipeline")
    snp.snpPipe(gtf, vcf, gos, snpAnnotOut)


# ARG PARSER
#-------------------------------------------------------------------------------
def parseArgs():

    parser = argparse.ArgumentParser(description='Unified Bioinformatician', add_help=False)
    parser.add_argument("-v","--verbose", action="store_true")
    subparsers = parser.add_subparsers(dest="cmd")

    snp_p = subparsers.add_parser("snps")
    snp_p.add_argument("-s","--vcf",required=True, help="Path to the vcf file with snps")
    snp_p.add_argument("-a","--gtf",required=True, help="Path to the gtf annotation file")
    snp_p.add_argument("-g","--gos",required=True, help="Path to the csv file with GO terms annotations (from Ensembl)")
    snp_p.add_argument("-o","--out",default="annotated_SNPs.csv", help="Path to the output file where annotated snps will be stored")
    snp_p.add_argument("-q","--gq_thres",default=20, help="The genotype quality threshold value to apply for genotyping")

    other_p = subparsers.add_parser("other")
    other_p.add_argument("name")

    args = parser.parse_args()

    return args

# MAIN
#-------------------------------------------------------------------------------

args = parseArgs()

if args.verbose:
     log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
else:
     log.basicConfig(format="%(levelname)s: %(message)s")

# The dictionnary of commands to launch depending on the "cmd"
# specified in the parser
cmds_d = dict(snps = snpAnalyse(args.gtf, args.vcf, args.gos, args.out))

# The command launch
cmds_d[args.cmd]
