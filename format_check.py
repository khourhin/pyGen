#-------------------------------------------------------------------------------
def isFasta(infile):
    with open(infile) as f:
        first_line = f.readline()

        if line.startswith(">"):
            return True
        else:
            return False

#-------------------------------------------------------------------------------
