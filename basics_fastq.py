import logging as log
import numpy
import basics_nuc_seq as bns


def get_casava_vers(fastq):
    """
    Check for the casava illumina pipeline version
    This only checks the first id line
    """

    with open(fastq, "r") as f:
        line1 = f.readline().rstrip()

        if line1.endswith("/1") or line1.endswith("/2"):
            return "<1.8"

        else:
            line1 = line1.split()[1]
            if line1.startswith("1") or line1.startswith("2"):
                return ">1.8"

            else:
                return "UNKNOWN"


def get_reads(fastq):
    """
    Return an iterator with all the reads
    """

    with open(fastq, "r") as f:

        try:
            for line in f:
                seq_id = line.strip()
                seq = f.next().strip()
                f.next()
                qual = f.next().strip()

                yield(seq_id, seq, qual)

        except StopIteration:
            raise IOError(
                "Can not read the fastq file: {}\nIs it a properly formatted fastq ?".format(fastq))


def fastq_stats(fastq):
    """
    Return stats about the fastq file
    """
    log.info("Loading data from: %s" % fastq)

    rds = get_reads(fastq)

    rds_len = []
    GC_content = []
    N_count = 0

    log.info("Computing statistics ...")

    for seq_id, seq, qual in rds:

        rds_len.append(len(seq))
        GC_content.append(bns.get_seq_GC(seq))
        N_count += seq.count("N")

    print("Reads statistics:")
    print("Total #reads: {}M".format(len(rds_len) / 1.0e6))
    print("Total #bases: {}G".format(sum(rds_len) / 1.0e9))
    print("Total #Ns: {}".format(N_count))
    print("Min length: {}".format(min(rds_len)))
    print("Max length: {}".format(max(rds_len)))
    print("Mean length: {}".format(numpy.mean(rds_len)))
    print("Standard deviation: {}".format(numpy.std(rds_len)))
    print("Mean GC content: {}".format(numpy.mean(GC_content)))


def filter_qual(fastq, fastq_out):
    """
    Perform the filter from Zhao 2011
    
    Filt 1:
    "reads that did not contain at least 41 Q20 bases among the
    first 51 cycles were removed."

    Filt 2:
    "Low quality (<Q20) 3prim" trimmed off 
    """
    # Quality dictionnaries:
    #---------------------------------------------------------------------------
    p33 = [ chr(i) for i in range(33,74) ]
    p33_d = { k:v for (k,v) in zip( p33, range(0,41) ) }

    reads = get_reads(fastq)

    with open(fastq_out, "w") as fout:
        for seq_id, seq, qual in reads:
            
            # translate phred+33 to num
            quals = [ p33_d[q] for q in qual ]

            # Filt 1:
            if len([ q for q in quals[0:52] if q >= 20 ]) < 41:
                pass

            else:

                #Filt 2:
                last_qual = quals[-1]
                while seq and last_qual < 20:
                    seq = seq[:-1]
                    qual = qual[:-1]
                    quals = quals[:-1]
                    last_qual = quals[-1]

                fout.write("{0}\n{1}\n{2}\n{3}\n".format(seq_id, seq, "+", qual))
                
#-------------------------------------------------------------------------------
def get_fq_ids(fastq):
    """
    Get all sequences id from a fastq file and return then along with the casava
    pipeline version
    """

    cas_vers = get_casava_vers(fastq)
    print "%s format: %s" % (fastq, cas_vers)
    
    id_lst = []
    with open(fastq, "r") as f:
        for line in f:
            # get ids
            id_lst.append( line.strip() )
            # pass other lines
            for i in range(3): next(f)
            
    return id_lst, cas_vers

#-------------------------------------------------------------------------------
def check_ID_first_field(fastq):
    """
    For testing, check if the first field (before first space) of the
    ID line can be used as an unique identifier in any case. (This 
    would make part of casava check obsolete I think ...)
    """

    id_list = []
    with open(fastq, "r") as f:
        for line in f:
            id_list.append(line.split()[0])
            for i in range(3): next(f)

    if len(id_list) == len(set(id_list)):
        print "First ID field seems unique"
    else:
        print "First ID field NOT unique"

#-------------------------------------------------------------------------------
def gen_reads(seq_d, rds_len, rds_n):
    """
    Create a number rds_n of paired reads of length rds_len from a seq
    dictionnary.
    """
    insert_size = 200

    for i in range(1,rds_n):
        seq_id = random.choice(seq_d.keys())
        # Define the maximum start point for the fragment to start
        max_start = len(seq_d[seq_id]) - rds_len * 2 - insert_size
        frag_start = random.randrange(0, max_start)
        frag_end = frag_start + rds_len * 2 + insert_size
        read_l = seq_d[seq_id][ frag_start : (frag_start + rds_len) ]
        read_r = seq_d[seq_id][ frag_end - rds_len: frag_end ]
        
        with open("my_reads_l.fas", "a") as f:
            f.write("@{}_{}/1\n{}\n+{}_{}/1\n{}\n".format(seq_id, i, read_l, seq_id ,i, "A"* rds_len))
        with open("my_reads_r.fas", "a") as f:
            f.write("@{}_{}/1\n{}\n+{}_{}/1\n{}\n".format(seq_id, i, read_r, seq_id ,i, "A"* rds_len))
