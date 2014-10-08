#! /usr/bin/python

#-------------------------------------------------------------------------------
def fasta_2_dict(fas_file, simple_ID=True):
    """
    Get a dictionnary from fasta file with:
    key = seq_id
    value = seq

    If simple_ID is True, then the seq_id in the dictionnary will be
    the original one truncated after the first white space.
    
    WORKING
    """
    with open(fas_file, "r") as f:
        seq_dict = {}
                
        # split each time a ">" in encountered        
        fasta = f.read().split('>')[1:]

        # split the strings in 3 items tuple: seq id, the sep, and the seq        
        tmp = [ seq.partition("\n") for seq in fasta ]

        # build a dictionnary with key=seq id, value=seq        
        for i in range( len( tmp ) ):
            
            if simple_ID:
                # for ex: for trinity output, split()[0] removes the len, path infos
                seq_dict[ tmp[i][0].split()[0] ] = tmp[i][2].replace( "\n", "").upper()
            else:
                seq_dict[ tmp[i][0] ] = tmp[i][2].replace( "\n", "").upper()

    return seq_dict

#-------------------------------------------------------------------------------
def get_these_seqs(seq_dict, seq_id_lst):
    """
    Yield tuples (seq_id, seq) from a LIST of selected ids
    """
    for seq_id in seq_id_lst:
            yield (seq_id, seq_dict[seq_id])

#-------------------------------------------------------------------------------
def write_to_fas(seq_dict, seq_id_lst, fas_out):
    """
    Write to the fas_out fasta file the selected seqs from the seq_ids LIST.
    Each seq line got a maximum of 80 characters
    """
    with open(fas_out, "w") as fout:

        for seq_id, seq in get_these_seqs(seq_dict, seq_id_lst):
            fout.write(">" + seq_id + "\n")

            # To cut seq lines each 80 characters
            while len(seq) > 80:
                fout.write(seq[:80] + "\n")
                seq =  seq[80:]
            fout.write(seq + "\n")

#-------------------------------------------------------------------------------
def print_summary(seq_dict):
    """
    Print out a summary of nucleotide sequences statistics
    """
    import basics_nuc_seq as bns
    
    all_GCs = bns.get_all_GCs(seq_dict)
    all_lens = bns.get_all_lens(seq_dict)
    n_seq = len(seq_dict)
    N50 = bns.get_N50(seq_dict)

    print "Total seqs:\t{0}".format(n_seq)
    print "Mean GC:\t{0}".format(sum(all_GCs.values()) / float(n_seq))
    print "Mean length:\t{0}".format(sum(all_lens.values()) / float(n_seq))
    print "Total length (in nucl):\t{0}".format(sum(all_lens.values()))
    print "N50:\t{0}".format(N50)
    print "Contigs in N50:\t{0}".format(len([x for x in all_lens.values() if x >= N50 ]))

#-------------------------------------------------------------------------------
if __name__ == "__main__":
    """
    FOR DEBUGGING
    """
    
    seq_dict = fasta_2_dict("demo_data/seq1.fas")
#    write_to_fas(seq_dict, ["seq1", "seq3"], "demo_data/myout")
    print_summary(seq_dict)
