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

if __name__ == "__main__":
    """
    FOR DEBUGGING
    """
    
    print fasta_2_dict("demo_data/seq1.fas", False)
