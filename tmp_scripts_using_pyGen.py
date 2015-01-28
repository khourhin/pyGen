#! /usr/bin/env python

import basics_fasta as bfas
import sys

# Filter sequences with a minimum length
fas_f = sys.argv[1]
min_len = sys.argv[2]

seqs_d = bfas.fasta_2_dict(fas_f)
seqs_i = bfas.filter_by_len(seqs_d, min_len)
bfas.write_as_fas(seqs_i)

