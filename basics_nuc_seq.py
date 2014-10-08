#! /usr/bin/python

#-------------------------------------------------------------------------------
def get_seq_GC(seq):
    """
    From a nucleotide sequence return the GC content
    """
    seq = seq.lower()
    n_C = seq.count("c")
    n_G = seq.count("g")
    
    GC_content = ( n_C + n_G ) / float(len(seq))
    return GC_content

#-------------------------------------------------------------------------------
def get_all_GCs(seq_dict):
    """
    From a seq_dict (fas2dict)
    Return all GCs as dictionnary with seqid as keys
    """
    all_GCs = {k: get_seq_GC( seq_dict[k] ) for k in seq_dict }
    
    return all_GCs

#-------------------------------------------------------------------------------
def get_all_lens(seq_dict):
    """
    From a seq_dict (fas2dict):
    Return all lengths as dictionnary with seqid as keys
    """
    all_lens = {k: len( seq_dict[k] ) for k in seq_dict }

    return all_lens

#-------------------------------------------------------------------------------
def get_N50(seq_dict):
    """
    From a seq_dict (fas2dict):
    Return the N50
    """
    import numpy 

    all_lens = get_all_lens(seq_dict)
    l_N50 = []

    for i in all_lens.values():
        for j in range(i):
            l_N50.append(i)

    N50 = numpy.median(l_N50)
    return N50

#-------------------------------------------------------------------------------
