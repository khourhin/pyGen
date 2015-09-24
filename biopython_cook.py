#! /usr/bin/env python
"""
A module with classical bioinfo tasks using Biopython.
So far:
- Getting a fasta file from a list of GIs
"""

from Bio import Entrez
from Bio import SeqIO
import csv
import sys

Entrez.email = 'ekornobis@gmail.com'

#-------------------------------------------------------------------------------
def getGIs(gis_csv, col_num=2):
    """
    Get the GIs from a csv file, specifying the number of the column
    where to find them (col index starting at 0)
    """
    gis = []

    with open(gis_csv, 'r') as f:
        next(f)
        for row in csv.reader(f):
            gis.append(row[col_num])
    return gis

#-------------------------------------------------------------------------------
def getFastaFromGIs(gis_list):
    """
    From a list of GIs (as string) get a multifasta
    """
    # Get data from NCBI
    handle = Entrez.efetch(db='protein', id=gis_list,
                           rettype='fasta', retmode='text')

    for record in SeqIO.parse(handle, 'fasta'):
        print '>' + record.id, record.description
        print record.seq

    handle.close()

#-------------------------------------------------------------------------------
if __name__ == '__main__':

    gis = getGIs(sys.argv[1])
    getFastaFromGIs(gis)
