#! /usr/bin/env python

import sys
import basics_fasta as bf

def filterFastaByIds(fasta, idsFile):
    seqs = bf.fasta_2_dict(fasta)
    ids_to_keep = []
    
    with open(idsFile, 'r') as f:
        ids_to_keep = [ seqid.strip() for seqid in f ]

    filt_seqs = bf.filter_by_id(seqs, ids_to_keep, False)

    bf.write_as_fas(filt_seqs)

def summarize_fasta(fas):
    seqs = bf.fasta_2_dict(fas)
    bf.make_summary(seqs)
    
    
#-------------------------------------------------------------------------------
if __name__ == '__main__':

    #filterFastaByIds(sys.argv[1], sys.argv[2])
    summarize_fasta( sys.argv[1] )





