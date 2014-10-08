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
        seq_d = {}
                
        # split each time a ">" in encountered        
        fasta = f.read().split('>')[1:]

        # split the strings in 3 items tuple: seq id, the sep, and the seq        
        tmp = [ seq.partition("\n") for seq in fasta ]

        # build a dictionnary with key=seq id, value=seq        
        for i in range( len( tmp ) ):
            
            if simple_ID:
                # for ex: for trinity output, split()[0] removes the len, path infos
                seq_d[ tmp[i][0].split()[0] ] = tmp[i][2].replace( "\n", "").upper()
            else:
                seq_d[ tmp[i][0] ] = tmp[i][2].replace( "\n", "").upper()

    return seq_d

#-------------------------------------------------------------------------------
def get_these_seqs(seq_d, seq_id_lst):
    """
    Yield tuples (seq_id, seq) from a LIST of selected ids
    """
    for seq_id in seq_id_lst:
            yield (seq_id, seq_d[seq_id])

#-------------------------------------------------------------------------------
def write_to_fas(seq_d, seq_id_lst, fas_out):
    """
    Write to the fas_out fasta file the selected seqs from the seq_ids LIST.
    Each seq line got a maximum of 80 characters
    """
    with open(fas_out, "w") as fout:

        for seq_id, seq in get_these_seqs(seq_d, seq_id_lst):
            fout.write(">" + seq_id + "\n")

            # To cut seq lines each 80 characters
            while len(seq) > 80:
                fout.write(seq[:80] + "\n")
                seq =  seq[80:]
            fout.write(seq + "\n")

#-------------------------------------------------------------------------------
def make_summary(seq_d, graphs_path):
    """
    Print out a summary of nucleotide sequences statistics
    """
    
    GCs_d = bns.get_all_GCs(seq_d)
    lens_d = bns.get_all_lens(seq_d)
    n_seq = len(seq_d)
    N50 = bns.get_N50(seq_d)

    # Stdout:
    print "Total seqs:\t{0}".format(n_seq)
    print "Total length (in nucl):\t{0}".format(sum(lens_d.values()))
    print "Mean GC:\t{0:.2f}".format(sum(GCs_d.values()) / float(n_seq))
    print "Mean length:\t{0:.2f}".format(sum(lens_d.values()) / float(n_seq))
    print "Standard deviation:\t{0:.2f}".format(numpy.std(lens_d.values()))
    print "Median length:\t{0}".format(numpy.median(lens_d.values()))
    print "Min length:\t{0}".format(min(lens_d.values()))
    print "Max length:\t{0}".format(max(lens_d.values()))
    print "N50:\t{0}".format(N50)
    print "Contigs in N50:\t{0}".format(len([x for x in lens_d.values() if x >= N50 ]))

    # Graphs:
    biog.plot_hist(lens_d, graphs_path + "Contig_lengths_histo.png" )
    
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    """
    FOR DEBUGGING
    """

    import basics_nuc_seq as bns
    import biographs as biog
    import numpy
    
    seq_d = fasta_2_dict("demo_data/AP1.fas")
    #    write_to_fas(seq_d, ["seq1", "seq3"], "demo_data/myout")
    make_summary(seq_d, "./")
    