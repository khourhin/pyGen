#-------------------------------------------------------------------------------
def get_lst_from_file(infile):
    """ 
    From a file with one entry per line, get a python list of these entries
    """

    with open(infile, "r") as f:
        elmt_l =  [ line for line in f]
        
    return elmt_l

#-------------------------------------------------------------------------------
