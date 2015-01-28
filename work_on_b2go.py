#! /usr/bin/python

import sys

#-------------------------------------------------------------------------------
def parse_b2go(annot_f, annot_d=None):
    """
    Parse the annotation file out of b2go i.e ".annot"
    Return a seqids dict of dicts with:
    GOs, ECs and annot as keys
    """
    with open(annot_f, "r") as f:

        for line in f:
            line = line.strip().split("\t")
            seq_id = line[0]

            # Create a dict entry if new seqid
            # This is weird (should it be created out of the for loop ?)
            if not annot_d:
                annot_d = {}

            if seq_id not in annot_d:
                annot_d[seq_id] = {}

            # Check and sort by annotation type
            for i in line[1:]:

                if i.startswith("GO:"):
                    # Method for creating or appending to existing list
                    annot_d[ seq_id ].setdefault("GOs", []).append(i)

                elif i.startswith("EC:"):
                    annot_d[ seq_id ].setdefault("ECs", []).append(i)

                    # Should not have more than 1 annot but for doublec hecking
                else:
                    annot_d[ seq_id ].setdefault("annot",[]).append(i)

    return annot_d

#-------------------------------------------------------------------------------
def parse_KEGG(kaas_f, annot_d):
    """
    Parse the output of the KEGG annotation from the KAAS webserver
    Return a annot dict
    """

    with open(kaas_f,"r") as f:

        for line in f:
            line = line.strip().split("\t")
            
            # If there is a kegg annotation (I think KEGGS should be unique)
            if not len(line) == 1:
                seqid = line[0]
                kegg = line[1]
                
                annot_d[seqid]["KEGG"] = kegg
                
    return annot_d

#-------------------------------------------------------------------------------
def print_annot(annot_d):
    """
    Print to std out:
    Tab delimited annotations:
    "SeqId","annot", "ECs", "GOs"
    """

    print "\t".join(["SeqId","annot","KEGG", "ECs", "GOs"])
    for seqid in annot_d:
        seq_d = annot_d[seqid]

        GOs = ",".join(seq_d["GOs"])
        annot = ",".join(seq_d["annot"])

        KEGG = seq_d["KEGG"] if "KEGG" in seq_d else "NA"
        ECs = ",".join(seq_d["ECs"]) if "ECs" in seq_d else "NA"

        print "\t".join([seqid, annot, KEGG, ECs, GOs])

#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------
if __name__ == "__main__":

    # Create an empty annotation dict
    annot_d = {}

    annot_d = parse_b2go(sys.argv[1], annot_d)
#    annot_d = parse_KEGG(sys.argv[2], annot_d)
    print_annot(annot_d)
    
    
