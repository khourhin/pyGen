import csv
import re


class Gene(object):
    """
    A gene representation for altanalyze
    """
    def __init__(self, row_dict):
        for k in row_dict:
            setattr(self, k, row_dict[k])
            self.junctions = {}

    def add_junction(self, junction):
        if junction.UID in self.junctions:
            print(junction.UID)
            print(self.Ensembl_gene)
            raise IOError('Duplicate Junction ids')
        else:
            self.junctions[junction.UID] = junction


class Junction(object):
    """
    A junction between two exons supported by read counts
    """
    def __init__(self, row_dict):
        for k in row_dict:
            setattr(self, k, row_dict[k])

        self.UID = tuple(row['UID'].split(':', 1)[1].split('-'))



    
     
   
count_file = "/home/ekornobis/analysis/muchardt/altAnalyze_visualization/Altanalyze_test/AltExpression/RNASeq/Hs/Hs_RNASeq_BS69_1_vs_Cont.ExpCutoff-5.0_average.txt"
annot_file = "/home/ekornobis/analysis/muchardt/altAnalyze_visualization/Altanalyze_test/ExpressionOutput/DATASET-AltAnalyze_Test2.txt"

study = {}

with open(annot_file) as infile:
    reader = csv.DictReader(infile, delimiter='\t')
    
    for row in reader:
        gene_id = row['Ensembl_gene']
        if gene_id in study:
            raise IOError('Duplicate Gene ids')
        else:
            study[gene_id] = Gene(row)

with open(count_file) as infile:
    reader = csv.DictReader(infile, delimiter='\t')

#### TO IMPROVE ####
# So far only use the lines with junction i.e matching "-" in first
# column this is a light regexp !!!
    
    for row in reader:
        gene_id = row['UID'].split(':')[0]
        gene = study[gene_id]

        if re.search("-", row['UID']):
            junc = Junction(row)
            gene.add_junction(junc)
