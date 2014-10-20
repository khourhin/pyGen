import subprocess

#-------------------------------------------------------------------------------
def format_db(fasta):
    """
    Format a fasta file for having them as db for blast+
    """

    cmd = [PATH_BLAST + 'makeblastdb', '-dbtype', 'nucl', '-in', fasta ]
    proc = subprocess.Popen( cmd, 
                             stdout=subprocess.PIPE)
    out, err = proc.communicate()

    return out, err

#-------------------------------------------------------------------------------
def do_blastN(query, db, outfile):
    """
    Launched a blast analysis for the [query] against the [db]
    """
    cmd = [PATH_BLAST + 'blastn', '-query', query, '-db', db, '-out', outfile,
           '-outfmt', '6', '-evalue', '1e-6' ]
    proc = subprocess.Popen( cmd, stdout=subprocess.PIPE)
    out, err = proc.communicate()

#-------------------------------------------------------------------------------

