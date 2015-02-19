#! /usr/bin/env python
# A library for xml blast results

#-------------------------------------------------------------------------------
def get_xml_res(blast_xml):
    """
    From a the blast_output.xml return a Blast Records object of results
    """

    from Bio.Blast import NCBIXML

    res = open(blast_xml, "r")
    blast_records = NCBIXML.parse(res)
    return blast_records
    res.close()

#-------------------------------------------------------------------------------
def filter_res(blast_rec, evalue):
    """
    Filter blast record
    """
    pass

#-------------------------------------------------------------------------------
def get_best_records(blast_recs):
    """
    From https://sites.google.com/site/xzhou82/Home/computer/python/biopython
    """

    for record in blast_recs:
        for align in record.alignments:
            hsp = align.hsps[0]
            # Remove spaces from the query ID
            queryID = record.query.split()[0]
            dbID = align.accession + align.title
            id_percent = float(hsp.identities) / hsp.align_length *100
            yield (queryID, dbID, id_percent, hsp.align_length, hsp.expect, hsp.score)
            break
    
#-------------------------------------------------------------------------------
def print_blast_records(blast_recs):
    """
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
    
    for blast_rec in blast_recs:
        for alignment in blast_rec.alignments:
            for hsp in alignment.hsps:

                # cp from blastxml_to_tabular.py (galaxy)
                mismatch = hsp.match.count(" ") + hsp.match.count("+") \
                           - hsp.query.count("-") - hsp.sbjct.count("-")
                gap_open = len(hsp.query.replace('-', ' ').split()) -1  \
                           + len(hsp.sbjct.replace('-', ' ').split()) -1 
                
                values = [ blast_rec.query.replace(" ", "_"),
                           alignment.title.replace(" ", "_"),
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
def print_gff(blast_recs):
    """
    IN DVPT: NOT WORKING
    TO CHECK MORE PROPERLY DEF OF GFF
    Return a gff3 file from a blast record

    PBS:
    STRAND DOES NOT seems to work (returning None None)
    """
    for blast_rec in blast_recs:
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

   import argparse

   parser = argparse.ArgumentParser(description="")
   parser.add_argument("infile",
                       help="xml input file")
#   parser.add_argument("outfile",
#                       help="output file")
#   parser.add_argument("-t","--test",
#                       help="an option",
#                       action='store_true')
   
   args = parser.parse_args()

   # Return result to standart out
   #blast_rec = get_xml_res(args.infile)
   #print_blast_records(blast_rec)

   #Only best records:
   blast_rec = get_xml_res(args.infile)
   print 'Query\tHit description\t%identity\tE-value\tScore'
   for i in get_best_records(blast_rec):
       print "\t".join([str(x) for x in i])


   #Return a gff
   # blast_rec = get_xml_res(args.infile)
   # print_gff(blast_rec)