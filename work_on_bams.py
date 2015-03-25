#! /usr/bin/python

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
def getCoverage(bam):
    """
    Get coverage for each genes using bedtools
    """
    log.info("Getting coverage for: %s" % bam)

    command = [ "bedtools", "genomecov", "-ibam", bam ]
    proc = subprocess.Popen(command)
    proc.communicate()
    print "ah"
    

#-------------------------------------------------------------------------------
if __name__ == "__main__":

    #getLocus(sys.argv[1],1,1,100)
    getCoverage(sys.argv[1])
    
