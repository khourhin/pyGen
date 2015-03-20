#! /usr/bin/python

import subprocess
import logging as log

#-------------------------------------------------------------------------------
def getLocus(bam, chromo, start, end):
    """
    Susbet a bam file with on a locus selected (using samtools)
    """

    command = [ "samtools", "view", "-bh", chromo + start + "-" + end ]
    proc = subprocess.Popen(command)
    proc.communicate()

    
    
    
    
