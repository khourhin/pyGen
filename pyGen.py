#! /usr/bin/env python
import format_check as fc
import basics_fasta as bfas
import basics_fastq as bfq
import blasting as blst
import common
import argparse
import logging as log

################################################################################
# A LIBRARY FOR GENETIC AND NGS DATA ANALYSIS
# UNDER CONSTRUCTION 
# Starting with Github ....
################################################################################

#-------------------------------------------------------------------------------
def sample_fasta(fasta, nseq):
    seq_d = bfas.fasta_2_dict(fasta)
    seq_i = bfas.get_random_seqs(seq_d, nseq)
    bfas.write_as_fas(seq_i)

#-------------------------------------------------------------------------------
def subset_fasta(fasta, seqids_l):
    """
    Get a subset of the fasta file with only the sequence which ids are in 
    ids_file.
    """
    seq_d = bfas.fasta_2_dict(fasta)
    seq_i = bfas.filter_by_id(seq_d, seqids_l)
    bfas.write_as_fas(seq_i)

#-------------------------------------------------------------------------------
def fasta_stats(fasta, graphs_path=None):

    seq_d = bfas.fasta_2_dict(fasta)
    bfas.make_summary(seq_d, graphs_path)

#-------------------------------------------------------------------------------
def blastN(query, db, outfile):
    """
    Perform a classical blastN of the query against the db. 
    Query and db should be in fasta
    outfile is the path to the blast output
    """
    blst.format_db(db)
    blst.do_blastN(query, db, outfile)

#-------------------------------------------------------------------------------
def file_stats(infile):
    
    if fc.isFasta(infile):
        print infile + " looks like a fasta file"
        fasta_stats(infile)

    elif fc.isFastQ(infile):
        print infile + " looks like a fastQ file"
        print "Casava version: " + bfq.get_casava_vers(infile)
        bfq.fastq_stats(infile)
#        bfq.check_ID_first_field(infile)
#        fastq_stats(infile)
    else:
        raise IOError("Formats others than Fasta or FastQ are not yet supported")
        
#-------------------------------------------------------------------------------
if __name__ == "__main__":


    # So far set by default:
    BLAST_OUT = "blast_out.tab"
    

    parser = argparse.ArgumentParser(description='Unified Bioinformatician')
    parser.add_argument("infile", help="the input file")
    parser.add_argument("-r","--random", type=int,
                        help="Return a random sample of 'R' sequences  ")
    parser.add_argument("-s", "--subset", help="Path to file with seqids to keep")
    parser.add_argument("-bN","--blastNdb", help="A fasta file to blast against")
    parser.add_argument("-S","--stats", action="store_true", help="Get stats")
    parser.add_argument("-Z","--tempScript", action="store_true", help="Tmp stuff")
    parser.add_argument("-v","--verbose", action="store_true", help="More logging info displayed")
    args = parser.parse_args()

    if args.infile:
        if args.verbose:
            log.basicConfig(format="%(levelname)s: %(message)s", level=log.DEBUG)
            log.info("Verbose is ON")
        else:
            log.basicConfig(format="%(levelname)s: %(message)s")
        
        if args.stats:
            file_stats(args.infile)

        if args.random:
            sample_fasta(args.infile, args.random)

        if args.subset:
            subset_fasta(args.infile, common.get_lst_from_file(args.subset))

        if args.blastNdb:
            blastN(args.infile, args.blastNdb, BLAST_OUT)

        if args.tempScript:
            # Filter fastq like Zhao 2011
            bfq.filter_qual(args.infile, "filtered_out.fq")
