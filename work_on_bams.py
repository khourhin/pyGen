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
if __name__ == "__main__":

getLocus(sys.argv[1],1,1,100)
    
    
    
