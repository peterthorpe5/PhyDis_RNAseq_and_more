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
    for line in f_in:
        input_count = input_count + 1
        if test_line(line):
            line = line.split()
            gene = line[0]
            goi_set.add(gene.rstrip())
    info = "input GOI list was %d genes" % input_count
    logger.info(info)
    info = "non-redundant GOI list was %d genes" % len(goi_set)
    logger.info(info)
    return(goi_set)


def parse_orthofinder(orthofinder, GOI_set, outfile,
                      logger):
    """fucn take in the  GOI_set
    if the gene in the orthofinder output check if gene in cluster. """
    f_in = open(orthofinder, "r")
    f_out = open(outfile, "w")
    matrix_found_count = 0
    for line in f_in:
        if test_line(line):
            elements = line.split()
            for protein in elements:
                for GOI in GOI_set:
                    # i thought it was a good idea to name the species
                    # Homo_NP_001001331.1 get rid of that. 
                    if "_" in protein:
                        protein = protein.replace("Homo_", "")
                        protein = protein.replace("Phy_", "")
                        protein = protein.replace("Mus_", "")
                    GOI = GOI.strip()
                    protein = protein.strip()
                    if GOI == protein:
                        out_data = "%s\t%s" % (GOI, line)
                        f_out.write(out_data)
    f_out.close()

           
if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python interog....py -c orthofinder.clusters.txt
--goi list_of_gene_of_interest -o outfile

"""

parser = OptionParser(usage=usage)

parser.add_option("-c", dest="orthofinder",
                  default="Orthogroups.txt",
                  help="Orthogroups.txt")

parser.add_option("--goi", dest="goi",
                  default=None,
                  help="file with a list of gene names",
                  metavar="FILE")
                  
parser.add_option("-o", dest="outfile",
                  default=None,
                  help="out file")



(options, args) = parser.parse_args()

orthofinder = options.orthofinder
goi = options.goi
outfile = options.outfile

logfile = "GOI_counts.log" 
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('get_GOI.py: %s'
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
    GOI_set = parse_goi(goi, logger)
    # parse the matrix

    
    parse_orthofinder(orthofinder, GOI_set, outfile,
                      logger)
            
