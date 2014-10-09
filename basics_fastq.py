#-------------------------------------------------------------------------------
def get_casava_vers(fastq):
    """
    Check for the casava illumina pipeline version
    This only checks the first id line
    """
    
    with open(fastq, "r") as f:
        line1 = f.readline().rstrip()
        
        if line1.endswith("/1") or line1.endswith("/2"):
            return "before_casava1.8"
        else:
            line1 = line1.split()[1]
            if line1.startswith("1") or line1.startswith("2"):
                return "after_casava1.8"
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
