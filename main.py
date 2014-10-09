#! /usr/bin/python
import basics_fasta as bfas
import format_check as fc

#-------------------------------------------------------------------------------
def fasta_stats(fasta, graphs_path=None):

    seq_d = bfas.fasta_2_dict(fasta)
    bfas.make_summary(seq_d, graphs_path)

#-------------------------------------------------------------------------------
def file_stats(infile):
    
    if fc.isFasta(infile):
        print infile + " looks like a fasta file"
        fasta_stats(infile)
    else:
        raise IOError("Other formats than Fasta are not yet supported")
        
#-------------------------------------------------------------------------------

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Unified Bioinformatician')
    parser.add_argument("infile", help="the input file")
    args = parser.parse_args()

    if args.infile:
        file_stats(args.infile)
