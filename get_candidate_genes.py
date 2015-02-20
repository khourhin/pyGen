#! /usr/bin/python

#-------------------------------------------------------------------------------
def extract_GO_from_b2go(annot):
    """
    return a dictionnary from an b2go.annot file
    keys: sequence names
    values: corresponding GO terms list
    """

    go_dict = {}
    annot_dict = {}
    with open(annot, "r") as f:
        for line in f:
            line = line.split()
            seq_id = line[0]
            go = line[1]
            # For now, skip the "EC" annotation
            if go.startswith("EC:"):
                continue
    
            if seq_id in go_dict:
                go_dict[ seq_id ].append(go)
            else:
                go_dict[ seq_id ] = [ go ]

            if line[2:]: # if there is an annotation
                annot_dict[ seq_id ] = " ".join(line[2:])
                
    return go_dict, annot_dict

#-------------------------------------------------------------------------------
def get_all_annot(annot_files_lst):
    """
    Get annotations from a list of files and return a dictionnary with:
    keys: file name
    values: list of annotations
    """
    import os
    
    total_annot = {}
    
    for f in annot_files_lst:
        go_dict, annot_dict = extract_GO_from_b2go(f)
        fname = os.path.basename(f)
        total_annot[fname] = annot_dict.values()

    return total_annot

#-------------------------------------------------------------------------------
def filter_single_annots(total_annot):
    """
    Filter a total_annot dictionnary to keep only annotations which are appearing
    only once within a sample
    """
    from collections import Counter
    
    for n, annot in total_annot.items():
        single_annots = [a for a, c in Counter(annot).items() if c == 1 ]
        total_annot[n] = single_annots

    return total_annot
    
#-------------------------------------------------------------------------------
def get_common_annots(total_annot):
    """
    From a total_annot dictionnary (filtered or not)),
    check which annotations are common to all 
    """
    from collections import Counter

    # Group all annotations to a single list:
    # set(j) in case filter_single_annots has not been used
    annots = [ i for j in total_annot.values() for i in set(j) ]

    # Get annots common to all (with a count equal to number of different samples)
    common_a = [ a for a, c in Counter(annots).items() if c == len(total_annot)]

    return common_a

#-------------------------------------------------------------------------------
def print_common_genes_summary(annot_files_lst):
    """
    Print a report about the annotations common to all files in a annot_files_lst
    """

    total_annot = get_all_annot(annot_files_lst)
    total_annot = filter_single_annots(total_annot)
    common_a = get_common_annots(total_annot)

    for annot in common_a:
        print annot

    print "Total common annotations: %d" % len(common_a)

#-------------------------------------------------------------------------------
def get_seq_names_by_annot(annot_lst, annot_dict):
    """
    Return the sequences names corresponding to the annotations in annot_lst
    provided in the annot_dict with the correspondance between seq names and annot
    """

    # Filter annot dict to get only the sequences with the annots in annot_lst
    # As well reverse the dict to get annotations as keys
    
    rev_annot = {annot_dict[k]: k for k in annot_dict
                  if annot_dict[k] in annot_lst }

    return rev_annot

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
def get_common_seqs(annot_files_lst, fastas_lst, outfolder):

    """
    From a list of annotations files and of fasta files (in the SAME ORDER)
    get a folder with a files with the sequences corresponding to an annotation
    common to all samples.
    """
    import os

    total_annot = get_all_annot(annot_files_lst)
    total_annot = filter_single_annots(total_annot)
    common_a = get_common_annots(total_annot)

    # Get annotations for each files and the seqnames corresponding to the common
    # annotations
    
    for fin, fasta in zip(annot_files_lst, fastas_lst):
        sample = os.path.basename(fasta)
        sample = os.path.splitext(sample)[0]
        
        print fasta
        go_dict, annot_dict = extract_GO_from_b2go(fin)
        rev_annot = get_seq_names_by_annot(common_a, annot_dict)

        seq_dict = fasta_2_dict(fasta)

    # Print out in a file for each annotation
        for annot in rev_annot:
            fname = annot.replace(" ", "_")
            fname = fname.replace("-", "_")
            
            seq_n = rev_annot[annot]
            
            with open(outfolder + fname, "a") as fout:
                fout.write( ">" + seq_n + "_" + sample + "\n")
                fout.write( seq_dict[seq_n] + "\n")
                
#    print get_all_seq_names_by_annot(common_a, annot_files_lst)
    
    
#-------------------------------------------------------------------------------
if __name__ == "__main__":

    import sys
    print "EDIT THE PROGRAM BEFORE LAUNCHING"
#    sys.exit()
    
    annots = [ "/home/tiennou/Desktop/transit/beetles/AP1.annot",
               "/home/tiennou/Desktop/transit/beetles/BB1.annot",
               "/home/tiennou/Desktop/transit/beetles/NC2.annot",
               "/home/tiennou/Desktop/transit/beetles/OO.annot",
               "/home/tiennou/Desktop/transit/beetles/thor.annot",
               "/home/tiennou/Desktop/transit/beetles/TM2.annot"
    ]
    fastas = [ "/home/tiennou/Desktop/transit/beetles/AP1_clus95.fas",
               "/home/tiennou/Desktop/transit/beetles/BB1_clus95.fas",
               "/home/tiennou/Desktop/transit/beetles/NC2_clus95.fas",
               "/home/tiennou/Desktop/transit/beetles/OO_clus95.fas",
               "/home/tiennou/Desktop/transit/beetles/Thor_clus95.fas",
               "/home/tiennou/Desktop/transit/beetles/TM2_clus95.fas"
    ]
    outfolder = "test3/"
    
    get_common_seqs(annots, fastas, outfolder)
