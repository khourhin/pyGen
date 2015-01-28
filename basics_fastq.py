#-------------------------------------------------------------------------------
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
                raise IOError("Your fastq seems neither casava1.8 or older")

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
