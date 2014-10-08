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
        fas_dict = {}
                
        # split each time a ">" in encountered        
        fasta = f.read().split('>')[1:]

        # split the strings in 3 items tuple: seq id, the sep, and the seq        
        tmp = [ seq.partition("\n") for seq in fasta ]

        # build a dictionnary with key=seq id, value=seq        
        for i in range( len( tmp ) ):
            
            if simple_ID:
                # for ex: for trinity output, split()[0] removes the len, path infos
                fas_dict[ tmp[i][0].split()[0] ] = tmp[i][2].replace( "\n", "").upper()
            else:
                fas_dict[ tmp[i][0] ] = tmp[i][2].replace( "\n", "").upper()

    return fas_dict

#-------------------------------------------------------------------------------
def get_these_seqs(fas_dict, seq_id_lst):
    """
    Yield tuples (seq_id, seq) from a LIST of selected ids
    """
    for seq_id in seq_id_lst:
            yield (seq_id, fas_dict[seq_id])

#-------------------------------------------------------------------------------
def print_to_fas(fas_dict, seq_id_lst, fas_out):
    """
    Print to the fas_out fasta file the selected seqs from the seq_ids LIST.
    Each seq line got a maximum of 80 characters
    """
    with open(fas_out, "w") as fout:

        for seq_id, seq in get_these_seqs(fas_dict, seq_id_lst):
            fout.write(">" + seq_id + "\n")

            # To cut seq lines each 80 characters
            while len(seq) > 80:
                fout.write(seq[:80] + "\n")
                seq =  seq[80:]
            fout.write(seq + "\n")

#-------------------------------------------------------------------------------
if __name__ == "__main__":
    """
    FOR DEBUGGING
    """
    
    fas_dict = fasta_2_dict("demo_data/seq1.fas")
    print_to_fas(fas_dict, ["seq1", "seq3"], "demo_data/myout")

