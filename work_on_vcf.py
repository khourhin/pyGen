#! /usr/bin/env python

# LOOKS LIKE WORKING, TO DOUBLE CHECK
# THIS WOULD NOT WORK WITH INDELS >>> TO FILTER INDELS FIRST
import sqlite3 as lite
import argparse

#-------------------------------------------------------------------------------
def get_genotype(all_b, geno_field, GQi, GQ_thres):

    # In case no data on genotype
    if geno_field.startswith("./."):
        return "NA"
    else:
        geno_f = geno_field.split(":")
        GQ = int(geno_f[GQi])

        if GQ >= GQ_thres:
            # Assume the genotype is the first field in the genotype field
            genotype = geno_f[0].split("/")
            # Form the genotype using the 0/1 notation as index
            genotype = [ all_b[int(x)] for x in genotype]
            return "/".join(genotype)
        else:
            return "FILT"

#--------------------------------------------------------------------------------
def get_SNP(line, GQ_thres):
    """
    Parse a SNP line from a vcf file
    """
    line = line.strip().split()
    chro = line[0]
    pos =  line[1]
    ref =  line[3]

    # Join all alternatives bases into a single string
    alt = "".join(line[4].split(","))
    all_b = ref + alt

    # Get the position of GQ in the INFO field
    GQi = line[8].split(":").index("GQ")

    geno_fields = line[9:]
    genotypes = [ get_genotype(all_b, x, GQi, GQ_thres) for x in geno_fields ]

    return (chro, pos, ref, all_b, "\t".join(genotypes))

#-------------------------------------------------------------------------------
def get_all_SNPs(vcf, GQ_thres):
    """
    Create a iterator over all SNPs from a vcf file
    """
    with open(vcf, "r") as f:
        for line in f:

            # Get he sample names
            if line.strip().startswith("#CHROM"):
                header = ["CHROM", "POS", "REF", "ALT"]
                header += line.strip().split()[9:]
                print "\t".join(header)
            # Skip comments
            if not line.strip().startswith("#"):
                yield get_SNP(line, GQ_thres)

#-------------------------------------------------------------------------------
def getGeneId(snp, gtf_db):
    """
    Fetch the gene name for the SNPs 
    from the gtf_db created by SNPs_annotator.py
    """
    chro = snp[0]
    pos = snp[1]

    with lite.connect(gtf_db) as con:
        con.text_factory = str # to not have unicode formatting
        cur = con.execute("""SELECT chromo, start, end, attributes from GTF
        WHERE chromo=? AND start < ? AND end > ?""",(chro, pos, pos))
        data = cur.fetchone()
        if data:
            gene_id = data[3].split('"')[1]
            return str(data[0]), str(data[1]), str(data[2]), gene_id
        else:
            return "NA", "NA", "NA", "NA"
    
#-------------------------------------------------------------------------------
def print_genotypes(vcf, GQ_thres, gtf_db=None):
    for snp in get_all_SNPs(vcf, GQ_thres):

        if gtf_db:
            gAnnot = getGeneId(snp, gtf_db)
            snp = gAnnot + snp

        print "\t".join( snp )

# MAIN
#-------------------------------------------------------------------------------

if __name__ == "__main__":

    # Not working yet with gq = 0 (tiny error)
    # This is quite slow (TODO make better use of the database)

    parser = argparse.ArgumentParser(description="Get genotypes from a VCF")
    parser.add_argument("VCF", help="A gtf file out of cufflinks to create a db")
    parser.add_argument("-t", "--GQ_thres", type=int,
                        help="A threshold for Genotype quality to be kept e.g. 20")
    parser.add_argument("-d", "--gtf_db",
                        help="The path the a gtf db built by SNPs_annotator.py")
    args = parser.parse_args()


    if args.VCF and args.GQ_thres:
        if args.gtf_db:
            print_genotypes(args.VCF, args.GQ_thres, args.gtf_db)
        else:
            print_genotypes(args.VCF, args.GQ_thres)
    else:
        raise IOError("Input should have VCF file and GQ threshold")
