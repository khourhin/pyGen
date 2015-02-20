#! /usr/bin/python

## The get get_Keggs function could be generalized to any annotation
## (with several columns, right now just using one columm)

import sys

#-------------------------------------------------------------------------------
def get_annot(annot_f):
    """
    From a tab delim files with:
    $1: seqID
    $2 .... $n annotations (b2go, kegg ...)
    Return a dictionnary
    """
    annot_d = {}
    with open(annot_f, "r") as f:
        for line in f:
            line = line.strip().split("\t")

            if len(line) == 1: # if no annot
                continue
            
            seqid = line[0] 
            annot = line[1:]
            

            if seqid in annot_d:
                raise IOError("Seems to have duplicate seqids: " + seqid)
            else:
                annot_d[seqid] = annot

    return annot_d
            
#-------------------------------------------------------------------------------
def get_clusters(clusters_f):
    """
    From a file out of corset with the clusters:
    $1: seqid
    $2: cluster_id
    Return a dictionny with:
    key: cluster ID
    values: list of seqids in the cluster
    """ 

    clus_d = {}
    with open(clusters_f, "r") as f:
        for line in f:
            seqid = line.split()[0]
            clusid = line.split()[1]

            if clusid in clus_d:
                clus_d[clusid].append(seqid)
            else:
                clus_d[clusid] = [seqid]

    return clus_d

#-------------------------------------------------------------------------------
def print_annotated_clusters(annot_d, clus_d):
    """
    Print tab delimited:
    $1 Clusterid
    $2 Annotations
    """

    for clus in clus_d:
        # Get all annots for each sequences in the cluster
        try:
            annots = set([ "\t".join(annot_d[seqid]) for seqid in clus_d[clus] ])
            print clus + "\t" + " ".join(annots)
        except KeyError:
            print clus + "\t" + "NA"
        
      
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    
    # USAGE:
    # work_on_corset KAAS_annot_file corset_clusters.txt

    # This is working now also with annot files created with
    # ~/bin/biopylib/work_on_b2go.py
    # BUT TO DOUBLE CHECK THIS PART

    annot_f = sys.argv[1]
    clus_f = sys.argv[2]
    
    annot_d = get_annot(annot_f)
    clus_d = get_clusters(clus_f)
    print_annotated_clusters(annot_d, clus_d)
