#! /usr/bin/python
import basics_fasta as bfas

#-------------------------------------------------------------------------------
def fasta_stats(fasta, graphs_path=None):

    seq_d = bfas.fasta_2_dict(fasta)
    bfas.make_summary(seq_d, graphs_path)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Unified Bioinformatician')
    parser.add_argument("infile", help="the input file")
    args = parser.parse_args()

    if args.infile:
        fasta_stats(args.infile)
