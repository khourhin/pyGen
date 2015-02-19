#! /usr/bin/python

# A library for extracting taxonomic information from a blast output
# USAGE:
# ./blast2tax.py <BLAST_OUTPUT>

#STATUS: SEEMS TO WORK
# ISSUES: lots of entries doesnt get taxonomy but seems to be from the
# updating of the db (records disapearing from genbank)

#-------------------------------------------------------------------------------
def get_GIs(blast_out):
    """
    With a blast output with only the best hits obtained by work_on_blast2.py
    Return a dictionnary of seq ids:
    k: seq_id
    v: GI

    The blast_out file should be with a header:
    Query	Hit	Hit description	E-value

    """
    
    import re
    
    p = re.compile('gi\|(\d*)\|') # define pattern for GI extraction
    id_dict = {}

    with open(blast_out, 'r') as f:

        next(f) # skip the header
        for line in f:
            m = re.search(p, line)
            GI = m.group(1)
            id_dict[line.split()[0]] = GI
            
    print "DONE: GIs parsed from the blast output"
    return id_dict

#-------------------------------------------------------------------------------
def get_taxID(id_dict):
    """
    From a set of GIs, get the tax id from gi_taxid_prot.dmp.gz, download from
    NCBI.
    Return a dictionnary with :
    k: GI
    v: taxid
    """
    import progressbar

    # Prepare progress bar
    num_lines = sum(1 for line in open('db/gi_taxid_prot.dmp', "r"))
    bar = progressbar.ProgressBar(maxval=num_lines,
                                  widgets=[progressbar.Bar('=', '[', ']'),
                                ' ', progressbar.Percentage()]) 

    GIs = set(id_dict.values())

    taxid_dict = {}
    found = 0
    
    with open('db/gi_taxid_prot.dmp', 'r') as ftax:
        print "Retreiving TaxIDs for each GI:"

        for i, line in enumerate(ftax):
            bar.update(i)
            line = line.strip().split('\t')
            
            if str(line[0]) in GIs:
                taxid_dict[str(line[0])] = int(line[1])
                found += 1

    bar.finish()
    print "DONE"
    print "Total GIs: %d" % len(GIs)
    print "Total found: %d" % found

    return taxid_dict

#-------------------------------------------------------------------------------
def get_taxonomy(taxid_dict):
    """
    Fetch from Entrez the complete lineage of all the entries in taxid_dict 
    Return a dict with:
    k: GI
    v: complete lineage
    """
    
    from Bio import Entrez
    Entrez.email = "ekornobis@gmail.com"

    tax_dict = {}
    for gi in taxid_dict:
        # Retrieve the taxonomy xml entry corresponding to the taxid
        rec = Entrez.efetch(db="taxonomy", id=taxid_dict[gi], retmode="xml")
        tax = Entrez.read(rec)
        tax_dict[gi] = tax[0]["Lineage"]

    print "DONE: Complete taxonomy retrived for each GI"
    return tax_dict

#-------------------------------------------------------------------------------
def get_tax_local(taxid_dict, ranks):
    """
    Fetch complete lineage locally (BETTER for huge queries)
    Return a dict with:
    k: GI
    v: complete lineage
    """    
    from cogent.parse.ncbi_taxonomy import NcbiTaxonomyFromFiles

    tree = NcbiTaxonomyFromFiles(open('db/nodes.dmp'), open('db/names.dmp'))
    root = tree.Root

    def get_lineage(node, my_ranks):
        ranks_lookup = dict([(r,idx) for idx,r in enumerate(my_ranks)])
        lineage = [None] * len(my_ranks)
        curr = node
        while curr.Parent is not None:
            if curr.Rank in ranks_lookup:
                lineage[ranks_lookup[curr.Rank]] = curr.Name
            curr = curr.Parent
        return lineage

    tax_dict = {}
    for gi in taxid_dict:
        # Get lineage for each (gi i.e taxid)
        try:
            node = tree.ById[taxid_dict[gi]]
            tax_dict[gi] = get_lineage(node, ranks)
        except KeyError:
            print "Cannot Fetch taxonomy for GI: " + gi

    print "DONE: Complete taxonomy retrived for each GI"
    return tax_dict
            
#-------------------------------------------------------------------------------
def print_tax_summary(tax_dict, id_dict, ranks):
    
    with open("tax_report.txt", "w") as f:
        f.write("Query\tHit_GI\t")
        f.write("\t".join(ranks) + "\n")
        for k in id_dict:
            gi = id_dict[k]
            try:
                tax = tax_dict[gi]
                # Change None values for string for printing
                tax = ['None' if v is None else v for v in tax]
                f.write("%s\t%s\t" % (k, gi))
                f.write("\t".join(tax) + "\n")
            except KeyError:
                f.write("%s\t%s\tTOCHECK\n" % (k, gi) )
        
    
#-------------------------------------------------------------------------------
if __name__ == '__main__':
    
    import sys
    # Can be edited
    ranks =  ['superkingdom','phylum','class','order','family','genus','species']

    id_dict = get_GIs(sys.argv[1])
    taxid_dict = get_taxID(id_dict)
    tax_dict = get_tax_local(taxid_dict, ranks)
    print_tax_summary(tax_dict, id_dict, ranks)

