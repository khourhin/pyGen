#! /usr/bin/env python
# A library for xml blast results

#-------------------------------------------------------------------------------
def get_xml_res(blast_xml):
    """
    From a the blast_output.xml return a Blast Records object of results
    """

    from Bio.Blast import NCBIXML

    with open(blast_xml, "r") as res:

        for blast_record in NCBIXML.parse(res):
            yield blast_record

#-------------------------------------------------------------------------------
def get_records(blast_xml, onlyBest=False):
    """
    From https://sites.google.com/site/xzhou82/Home/computer/python/biopython
    """
    for record in get_xml_res(blast_xml):

        for align in record.alignments:

            queryID = record.query.split()[0]
            dbID = align.hit_def.replace(" ", "_")
            #old version (don't know exactly what align.accession is refering to)
            #dbID = align.accession + align.title
            
            if onlyBest:
                hsp = align.hsps[0]
                yield (queryID, dbID, hsp)
                
            else:
                for hsp in align.hsps:
                    yield (queryID, dbID, hsp)

#-------------------------------------------------------------------------------
def print_recs(hsp_i):
    """
    from a tuple hsp iterator of the form:
    (queryID, dbID, <hsp object>) # like in get_best_records

    Print a blast records in tab delimited format with fields as:

    ====== ========= ============================================
    Column NCBI name Description
    ------ --------- --------------------------------------------
         1 qseqid    Query Seq-id (ID of your sequence)
         2 sseqid    Subject Seq-id (ID of the database hit)
         3 pident    Percentage of identical matches
         4 length    Alignment length
         5 mismatch  Number of mismatches
         6 gapopen   Number of gap openings
         7 qstart    Start of alignment in query
         8 qend      End of alignment in query
         9 sstart    Start of alignment in subject (database hit)
        10 send      End of alignment in subject (database hit)
        11 evalue    Expectation value (E-value)
        12 bitscore  Bit score
    ====== ========= ============================================
    """

    for hsp in hsp_i:
        queryID = hsp[0]
        dbID = hsp[1]
        hsp = hsp[2]
        
        mismatch = hsp.match.count(" ") + hsp.match.count("+") \
                   - hsp.query.count("-") - hsp.sbjct.count("-")
        gap_open = len(hsp.query.replace('-', ' ').split()) -1  \
                   + len(hsp.sbjct.replace('-', ' ').split()) -1 
        
        values = [ queryID,
                   dbID,
                   "%.2f"% (float(hsp.identities) / hsp.align_length *100),
                   mismatch,
                   gap_open,
                   hsp.align_length,
                   hsp.query_start,
                   hsp.query_end,
                   hsp.sbjct_start,
                   hsp.sbjct_end,
                   hsp.expect,
                   hsp.score ]
        
        # Convert all valuess to strings
        values = [ str(x) for x in values ]
        
        print "\t".join(values)
        
    
#-------------------------------------------------------------------------------
def print_gff(blast_xml):
    """
    IN DVPT: NOT WORKING
    TO CHECK MORE PROPERLY DEF OF GFF
    Return a gff3 file from a blast record

    PBS:
    STRAND DOES NOT seems to work (returning None None)
    """
    for blast_rec in get_xml_res(blast_xml):
        
        for alignment in blast_rec.alignments:
            for hsp in alignment.hsps:

                attr = "ID:%s" % alignment.title.replace(" ", "_"),

                values = [ blast_rec.query.replace(" ", "_"),
                           "blast",
                           "transcript",
                           hsp.query_start,
                           hsp.query_end,
                           hsp.expect,
                           hsp.strand,
                           ".",
                           attr ]

                # Convert all valuess to strings
                values = [ str(x) for x in values ]

                print "\t".join(values)
                
#-------------------------------------------------------------------------------
# THE MAIN ---------------------------------------------------------------------
#-------------------------------------------------------------------------------

if __name__ == "__main__":

    import sys
    
    # Get best hits:
    xml = sys.argv[1]

    recs = get_records(xml, onlyBest=True)
    print_recs(recs)
