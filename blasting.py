#! /usr/bin/env python

import sys
import subprocess

# WARNING: Only best is not properly working !!!
# '-max_target_seqs' dows not report only one hit per query

PATH_BLAST = '/home/tiennou/Documents/Taff/softwares/RNA-seq/blast+/bin/'

#-------------------------------------------------------------------------------
def format_db(fasta, dbtype='nucl'):
    """
    Format a fasta file for having them as db for blast+
    dbtype can be 'nucl' or 'prot'
    """
    print "Formatting DB"
    cmd = [PATH_BLAST + 'makeblastdb', '-dbtype', dbtype, '-in', fasta ]
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE)
    out, err = proc.communicate()

    return out, err

#-------------------------------------------------------------------------------
def do_blastN(query, db, outfile, onlyBest=False, cpus=1, fmt=6):
    """
    Launched a blastn analysis for the [query] against the [db]
    """
    print "Starting BlastN"
    cmd = [PATH_BLAST + 'blastn', '-query', query, '-db', db, '-out', outfile,
           '-outfmt', str(fmt), '-evalue', '1e-6', '-num_threads', str(cpus) ]
    if onlyBest:
        cmd = cmd + [ '-max_target_seqs', '1' ]
    print " ".join(cmd)
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE)
    out, err = proc.communicate()

#-------------------------------------------------------------------------------
def do_blastP(query, db, outfile, onlyBest=False, cpus=1, fmt=6):
    """
    Launched a blastp analysis for the [query] against the [db]
    """
    print "Starting BlastP"
    cmd = [PATH_BLAST + 'blastp', '-query', query, '-db', db, '-out', outfile,
           '-outfmt', str(fmt), '-evalue', '1e-6', '-num_threads', str(cpus) ]
    if onlyBest:
        cmd = cmd + [ '-max_target_seqs', '1' ]
    print " ".join(cmd)
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE)
    out, err = proc.communicate()

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    query = sys.argv[1]
    db = sys.argv[2]
    outfile = sys.argv[3]

#    format_db(db, 'prot')
    do_blastP(query, db, outfile, True, 2)

#    format_db(db, 'nucl')
#    do_blastN(query, db, outfile, True, 2)
