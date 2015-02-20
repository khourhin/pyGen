#! /usr/bin/python

# SEEMS TO WORK

import sys
import subprocess
import glob
import os

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
def format_all_db(fastas):
    """
    Format all fasta files for blast+
    """

    for fasta in fastas:
        out, err = format_db(fasta)
        
#-------------------------------------------------------------------------------
def do_blast(query, db, outfile):
    """
    Launched a blast analysis for the [query] against the [db]
    """
    cmd = [PATH_BLAST + 'blastn', '-query', query, '-db', db, '-out', outfile,
           '-outfmt', '6', '-evalue', '1e-6' ]
    proc = subprocess.Popen( cmd, 
                             stdout=subprocess.PIPE)
    out, err = proc.communicate()

#-------------------------------------------------------------------------------
def do_all_blast(fastas, outDir):
    """
    Launch all the blast (one search and its reciproc)
    """
    
    # Do the blast searches with first fasta against all other fastas
    # as dbs (and the reciproq search)

    out_lst = []
    query = fastas[0]
    dbs = fastas[1:]
    
    for db in dbs:
        print "Blasting " + query + " VS " + db
        # Get a mix of the 2 filenames
        outfile = outDir + "_VS_".join(
            [ os.path.basename(x) for x in (query, db) ])
        do_blast(query, db, outfile)
        do_blast(db, query , outfile + "_reci")
        out_lst.append(outfile)
        
        # Get list of blast out files 
    return out_lst

#-------------------------------------------------------------------------------
def get_hits(blast_out):
    """
    Store the best hits blast results in a dictionnary 
    """

    import sys

    hit_dict = {}

    with open(blast_out, "r") as f:
        for line in f:
            res = line.split()

            if res[0] in hit_dict:
                if res[-1] > hit_dict[res[0]][1]:
                    hit_dict[res[0]] = (res[1], res[-1])
            else:
                hit_dict[res[0]] = (res[1], res[-1])

        hit_dict = {key: value[0] for (key, value) in hit_dict.items() }
        return hit_dict

#-------------------------------------------------------------------------------
def check_all_reciproq(outdir, out_lst):
    """
    Check if the reciproq blasts gives the corresponding hit and return a 
    dictionnary with all good reciprocal hits
    """

    total_dict = {}
    # for each blast outfile ...
    for b_out in out_lst:
    
        b1_dict = get_hits(b_out)
        b2_dict = get_hits(b_out + "_reci") # ... and their reciproque search

        # Check if hits are "reciprocal"
        for query in b1_dict:
            hit = b1_dict[query]
            rec_hit = b2_dict.get(hit) #get is for not raising exception if no hit
        
            if query == rec_hit:
                b_out = os.path.basename(b_out) 
                key_name = (b_out.split("_VS_")[0], query)
                val_name = (b_out.split("_VS_")[1], hit)
                
                if key_name in total_dict:
                    total_dict[key_name].append(val_name) 
                else:
                    total_dict[key_name] = [ val_name ]

    return total_dict

#-------------------------------------------------------------------------------
def make_summary(total_dict, out_lst, outfile):

    final_dict = {}
    with open(outfile, "w") as out:

        for k in total_dict:

            # Check for homologs present in all datasets
            if len(total_dict[k]) == len(out_lst):
                final_dict[k] = total_dict[k]
                out.write(k[0] + ":" + k[1])
                for i in final_dict[k]:
                    out.write(" " + ":".join(i))
                out.write("\n")

    return final_dict

#-------------------------------------------------------------------------------
def fasta_2_dict(fasta_file):
    """
    Get a dictionnary from fasta file with:
    key = seq_id
    value = seq
    
    WORKING
    """
    
    with open(fasta_file, "r") as f:
        
        seq_dict = {}
        
        fasta = f.read().split('>')[1:]
        # split each time a ">" in encountered
        
        tmp = [ seq.partition("\n") for seq in fasta ]
        # split the strings in 3 items tuple: seq id, the sep, and the seq
        
        for i in range( len( tmp ) ):
            seq_dict[ tmp[i][0].split()[0] ] = tmp[i][2].replace( "\n", "").upper()
            # build a dictionnary with key=seq id, value=seq
            # made for trinity output: split()[0] removes the len, path infos
            
    return seq_dict

#-------------------------------------------------------------------------------
def make_alignments(fastas, final_dict, fasta_dir):

    all_fas_dict = {}
    f_count = 1

    # Get all fastas file as dictionnaries in a dictionnary
    for fas in fastas:
        print fas
        all_fas_dict[os.path.basename(fas)] = fasta_2_dict(fas)

    # Write each "homologs" sequences to a file
    for key in final_dict:
        with open(fasta_dir + str(f_count) + ".fas", "w") as f_out:
            # Write the "intial query" sequence
            f_out.write(">" + key[1] + ":" + key[0] + "\n" )
            f_out.write( all_fas_dict[key[0]][key[1]] + "\n")

            # Write others hit sequences
            for f_in, seq_id in final_dict[key]:
                f_out.write(">" + seq_id + ":" + f_in + "\n")
                f_out.write( all_fas_dict[f_in][seq_id] + "\n")
        f_count += 1
            
#-------------------------------------------------------------------------------
if __name__ == '__main__':

    # EVALUE IS SET TO 1E-6
    PATH_BLAST = '/home/tiennou/Documents/Taff/softwares/RNA-seq/blast+/bin/'
    outfile = "summary_test.txt"

    # IO:
    output_dir = sys.argv[1] + "/"
    fas_out_dir = output_dir + 'fastas/'
    fastas = sys.argv[2:]

    # Create output folders
    os.makedirs(fas_out_dir)

    format_all_db(fastas)
    out_lst = do_all_blast(fastas, output_dir)
    total_dict = check_all_reciproq(output_dir, out_lst)
    final_dict = make_summary(total_dict, out_lst, outfile)

    make_alignments(fastas, final_dict, fas_out_dir)
    
