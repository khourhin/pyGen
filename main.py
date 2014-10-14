#! /usr/bin/python
import format_check as fc
import basics_fasta as bfas
import basics_fastq as bfq

#-------------------------------------------------------------------------------
def fasta_stats(fasta, graphs_path=None):

    seq_d = bfas.fasta_2_dict(fasta)
    bfas.make_summary(seq_d, graphs_path)

#-------------------------------------------------------------------------------
def file_stats(infile):
    
    if fc.isFasta(infile):
        print infile + " looks like a fasta file"
        fasta_stats(infile)

    elif fc.isFastQ(infile):
        print infile + " looks like a fastQ file"
        print "Casava version: " + bfq.get_casava_vers(infile)
        bfq.check_ID_first_field(infile)
#        fastq_stats(infile)
    else:
        raise IOError("Formats others than Fasta oro FastQ are not yet supported")
        
#-------------------------------------------------------------------------------

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Unified Bioinformatician')
    parser.add_argument("infile", help="the input file")
    args = parser.parse_args()

    if args.infile:
        file_stats(args.infile)
