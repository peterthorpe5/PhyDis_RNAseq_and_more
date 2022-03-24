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
from modules.database import parse_NCBI_gffcol9_info


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
            
