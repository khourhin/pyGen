#! /usr/bin/env python

import subprocess
import logging as log
import sys

#-------------------------------------------------------------------------------
def getLocus(bam, chromo, start, end):
    """
    Susbet a bam file with on a locus selected (using samtools)
    """

    command = [ "samtools", "view", "-bh", chromo + start + "-" + end ]
    proc = subprocess.Popen(command)
    proc.communicate()


#-------------------------------------------------------------------------------
def getCoverage(bam, bedout):
    """
    Get coverage for each genes using bedtools
    Create an histogram using bedtools and then print to stdout
    the average coverage for each genes/transcripts
    """
    log.info("Getting coverage for: %s" % bam)

    # Get coverage histogram from bedtools
    with open(bedout, "w") as fout:
    
        command = [ "bedtools", "genomecov", "-ibam", bam ]
        proc = subprocess.Popen(command, stdout=subprocess.PIPE)
        fout.write( proc.communicate()[0] )

    with open(bedout, "r") as fin:
        cov_d = {}
        for entry in fin:
            entry = entry.split("\t")
            ref_seq = entry[0]

            cov_d[ref_seq] = cov_d.setdefault(ref_seq, 0) \
                             + int(entry[1]) * int(entry[2]) / float(entry[3])

    for k in cov_d:
        print "%s\t%0.1f" % (k, cov_d[k])
        
#-------------------------------------------------------------------------------
if __name__ == "__main__":

    bedout = "coverage_from_bed.tab"
    
    #getLocus(sys.argv[1],1,1,100)
    getCoverage(sys.argv[1], bedout)
    
