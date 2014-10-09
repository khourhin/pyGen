#! /usr/bin/python

#-------------------------------------------------------------------------------
def fasta_stats(fasta, graphs_path=None)

    import basics_nuc_seq as bns
    import biographs as biog
    import numpy
    
    seq_d = fasta_2_dict("demo_data/AP1.fas")
    make_summary(seq_d, graphs_path)



#-------------------------------------------------------------------------------

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Unified Bioinformatician')
    parser.add_argument("infile", help="the input file")
    args = parser.parse_args()

    if args.infile:
        fasta_stats(args.infile)
