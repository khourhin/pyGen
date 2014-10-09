#-------------------------------------------------------------------------------
def isFasta(infile):
    with open(infile) as f:
        line1 = f.readline()

        if line1.startswith(">"):
            return True
        else:
            return False

#-------------------------------------------------------------------------------
def isFastQ(infile):
    with open(infile) as f:
        line1 = f.readline()

        if line1.startswith("@"):
            return True
        else:
            return False

#-------------------------------------------------------------------------------    
