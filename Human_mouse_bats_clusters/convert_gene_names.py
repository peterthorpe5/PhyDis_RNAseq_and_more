# script to interogate clusters for GOI - convert gene names though
# author Pete Thorpe 2022 march
# imports

import os
import sys
from optparse import OptionParser
import datetime
import logging
import logging.handlers
import time
from collections import defaultdict


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    if "row" in line: # this is line 1
        return False
    return line


def parse_goi(goi, logger):
    """func take in a file with a list og gene names. 
    returns a set"""
    f_in = open(goi, "r")
    input_count = 0
    goi_set = set([])
    goi_to_line_dict = dict()
    for line in f_in:
        input_count = input_count + 1
        if test_line(line):
            line_str = line
            line = line.split()
            gene = line[0]
            goi_set.add(gene.rstrip())
            goi_to_line_dict[gene] = line_str
    info = "input GOI list was %d genes" % input_count
    logger.info(info)
    info = "non-redundant GOI list was %d genes" % len(goi_set)
    logger.info(info)
    return(goi_set, goi_to_line_dict)


def parse_NCBI_gffcol9_info(infile):
    """func to get the required info from a preparpared info file, as such

    zcat GCF_000001635.27_GRCm39_genomic.gff.gz | grep "CDS"
    | cut -f9 | uniq > mus_gene_to_symblo.info
    
    This will result sin output like this:
    ID=cds-XP_006503333.1;Parent=rna-XM_006503270.5;Dbxref=GeneID:66260,
    Genbank:XP_006503333.1,MGI:MGI:1913510;Name=XP_006503333.1;gbkey=CDS;
    gene=Tmem54;product=transmembrane protein 54 isoform X1;
    protein_id=XP_006503333.1

    we are interested in saving the protein_id= cds- gene= Name= product=

    The AA file does not match the gene ids used by the group, so
    this info is required for a traslation.
    e.g. human "A3GALT2" results in this:
    
    ID=cds-NP_001073907.1;Parent=rna-NM_001080438.1;Dbxref=CCDS:CCDS60080.1,
    Ensembl:ENSP00000475261.1,GeneID:127550,Genbank:NP_001073907.1,
    HGNC:HGNC:30005;Name=NP_001073907.1;gbkey=CDS;gene=A3GALT2;
    product=alpha-1%2C3-galactosyltransferase 2;protein_id=NP_001073907.1;
    tag=MANE Select
    """
    # set up some default dicts
    gene_to_prot = defaultdict(str)
    prot_id_to_gene = defaultdict(int)
    gene_to_product = defaultdict(str)
    prot_id_to_product = defaultdict(int)
    
    f_in = open(infile, "r")
    input_count = 0
    goi_set = set([])
    for line in f_in:
        input_count = input_count + 1
        if test_line(line):
            try:
                cds_info = line.split("ID=cds-")[1]
                cds_info = cds_info.split(";")[0]
            except:
                # no cds info
                continue
            try:
                # no gene info
                gene = line.split("gene=")[1]
                gene = gene.split(";")[0]
            except:
                continue
            try:
                product = line.split("product=")[1]
                product = product.split(";")[0]
            except:
                # no product info
                product = ""
            # print("cds_info, gene, product")
            # print(cds_info, gene, product)
            gene_to_prot[gene] = cds_info
            prot_id_to_gene[cds_info] = gene
            gene_to_product[gene] = product
            prot_id_to_product[cds_info] = product
    return gene_to_prot, prot_id_to_gene, gene_to_product, prot_id_to_product



def convert_names(GOI_set, goi_to_line_dict,
                  outfile, gene_to_prot,
                  prot_id_to_gene, gene_to_product,
                  prot_id_to_product,
                  logger):
    """fucn take in the counts.matrix and GOI_set
    if the gene in the matrix is in the set,
    prints it to the file. """
    f_out = open(outfile, "w")
    for GOI in GOI_set:
        print(GOI)
        # conversion.
        prot_name = gene_to_prot[GOI]
        original_line = goi_to_line_dict[GOI]
        new_line = original_line.replace(GOI, prot_name)
        f_out.write(new_line)
    f_out.close()



if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python interog....py -c orthofinder.clusters.txt
--goi list_of_gene_of_interest -o outfile

"""

parser = OptionParser(usage=usage)


parser.add_option("--goi", dest="goi",
                  default="GOI.txt",
                  help="file with a list of gene names",
                  metavar="FILE")

parser.add_option("-i", dest="info",
                  default=None,
                  help="col9 from the gff file",
                  metavar="FILE")
                  
parser.add_option("-o", dest="outfile",
                  default=None,
                  help="out file.matrix")



(options, args) = parser.parse_args()

goi = options.goi
info = options.info
outfile = options.outfile

logfile = "convert_gene_names.log" 
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('convert_gene_names.py: %s'
                               % time.asctime())
    logger.setLevel(logging.DEBUG)
    err_handler = logging.StreamHandler(sys.stderr)
    err_formatter = logging.Formatter('%(levelname)s: %(message)s')
    err_handler.setFormatter(err_formatter)
    logger.addHandler(err_handler)
    try:
        logstream = open(logfile, 'w')
        err_handler_file = logging.StreamHandler(logstream)
        err_handler_file.setFormatter(err_formatter)
        # logfile is always verbose
        err_handler_file.setLevel(logging.INFO)
        logger.addHandler(err_handler_file)
    except:
        logger.error("Could not open %s for logging" %
                     logfile)
        sys.exit(1)
    # get gene of interests
    GOI_set, goi_to_line_dict = parse_goi(goi, logger)
    # parse the matrix
    gene_to_prot, prot_id_to_gene, gene_to_product, \
            prot_id_to_product = parse_NCBI_gffcol9_info(info)
    # lets have a look at the first few elements of the dicts:
##    for mydict in [gene_to_prot, prot_id_to_gene, \
##                   gene_to_product, prot_id_to_product]:
##        first2pairs = {k: mydict[k] for k in list(mydict)[:2]}
##        print(first2pairs)
##    print("... testing converting OR4F5 to prot name. should be  NP_001005484.2...\n\n")
##    test_prot = gene_to_prot["OR4F5"]
##    print(" result = ", test_prot)
    
    convert_names(GOI_set, goi_to_line_dict,
                  outfile, gene_to_prot,
                  prot_id_to_gene, gene_to_product,
                  prot_id_to_product,
                  logger)
            
