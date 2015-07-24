#! /usr/bin/env python

import sys
import subprocess

PATH_BLAST = '/home/tiennou/Documents/Taff/softwares/RNA-seq/blast+/bin/'

#-------------------------------------------------------------------------------
def format_db(fasta, dbtype='nucl'):
    """
    Format a fasta file for having them as db for blast+
    dbtype can be 'nucl' or 'prot'
    """

    cmd = [PATH_BLAST + 'makeblastdb', '-dbtype', dbtype, '-in', fasta ]
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE)
    out, err = proc.communicate()

    return out, err

#-------------------------------------------------------------------------------
def do_blastN(query, db, outfile, onlyBest=False):
    """
    Launched a blastn analysis for the [query] against the [db]
    """
    cmd = [PATH_BLAST + 'blastn', '-query', query, '-db', db, '-out', outfile,
           '-outfmt', '6', '-evalue', '1e-6' ]
    if onlyBest:
        cmd = cmd + [ '-max_target_seqs', '1' ]
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE)
    out, err = proc.communicate()

#-------------------------------------------------------------------------------
def do_blastP(query, db, outfile, onlyBest=False):
    """
    Launched a blastp analysis for the [query] against the [db]
    """
    cmd = [PATH_BLAST + 'blastp', '-query', query, '-db', db, '-out', outfile,
           '-outfmt', '6', '-evalue', '1e-6']
    if onlyBest:
        cmd = cmd + [ '-max_target_seqs', '1' ]
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE)
    out, err = proc.communicate()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    query = sys.argv[1]
    db = sys.argv[2]
    outfile = sys.argv[3]

#    format_db(db, 'prot')
    do_blastP(query, db, outfile, True)
