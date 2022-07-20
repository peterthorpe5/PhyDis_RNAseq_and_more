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
from modules.database import parse_NCBI_gffcol9_info, test_line




def parse_file(infile, gene_to_up_or_down, gene_to_LFC, logger):
    """func to parse the de results and writes the sig results"""
    count = 0
    name_set = set([])
    data = open (infile, "r")
    gene_names = set([])

    for line in data:
        if test_line(line):
            test = line.split()
            if len(test) !=3 : continue
            name, logfc, FDR =  line.split()
            logfc = float(logfc)
            FDR = float(FDR)
            if name in name_set:
                # print("duplicate", name)
                continue
            # add the name to the set        
            name_set.add(name)
            if logfc > 0:
                gene_to_up_or_down[name] = "up"
                gene_to_LFC[name] = logfc
            # convert the negative to positive for easy testing
            if logfc < 0:
                gene_to_up_or_down[name] = "down"
                gene_to_LFC[name] = logfc
    return gene_to_up_or_down, gene_to_LFC


def parse_orthofinder(orthofinder, outfile, hu_gene_to_prot,
                      mo_gene_to_prot,
                      gene_to_up_or_down, gene_to_LFC,
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
                # i thought it was a good idea to name the species
                # Homo_NP_001001331.1 get rid of that. 
                if "_" in protein:
                    protein = protein.replace("Homo_", "")
                    protein = protein.replace("Phy_", "")
                    protein = protein.replace("Mus_", "")
                protein = protein.strip()
                LFC = gene_to_LFC[protein]
                print(protein)
                print(LFC)
                if float(LFC) > 0.4:
                    out_data = "%s (LOGFC = %s) " % (protein, LFC)
                    line = line.replace(protein, out_data)
            f_out.write(line)
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
                 

parser.add_option("--in1", dest="in1",
                  default="/storage/home/users/sonia_vernes/human_mouse_bat_clustering/filtered_by_logFC_0.40_FDR_0.01/human_prot_id_LOGFC_0.40_FDR_0.01",
                  help="human filtered RNAseq results")

parser.add_option("--in2", dest="in2",
                  default="/storage/home/users/sonia_vernes/human_mouse_bat_clustering/filtered_by_logFC_0.40_FDR_0.01/mouse_prot_id_LOGFC_0.40_FDR_0.01",
                  help="mouse filtered RNAseq results")

parser.add_option("--in3", dest="in3",
                  default="/storage/home/users/sonia_vernes/human_mouse_bat_clustering/filtered_by_logFC_0.40_FDR_0.01/phydis_LOGFC_0.40_FDR_0.01",
                  help="phydis filtered RNAseq results")

parser.add_option("--hu", dest="human",
                  default="human_gene_to_symblo.info",
                  help="human_gene_to_symblo.info")

parser.add_option("--mu", dest="mouse",
                  default="mus_gene_to_symblo.info",
                  help="mus_gene_to_symblo.info")                  
parser.add_option("-o", dest="outfile",
                  default="rewite_test",
                  help="out file")



(options, args) = parser.parse_args()

in1 = options.in1
in2 = options.in2
in3 = options.in3

human = options.human
mouse = options.mouse
orthofinder = options.orthofinder
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
    # parse the matrix
    hu_gene_to_prot, hu_prot_id_to_gene, hu_gene_to_product, \
            hu_prot_id_to_product = parse_NCBI_gffcol9_info(human)

    mo_gene_to_prot, mo_prot_id_to_gene, mo_gene_to_product, \
            mo_prot_id_to_product = parse_NCBI_gffcol9_info(mouse)

    # set up some dictionaries to capture the data. 
    gene_to_up_or_down = defaultdict(str)
    gene_to_LFC = defaultdict(float)
    # add the human data
    gene_to_up_or_down, gene_to_LFC = parse_file(in1, gene_to_up_or_down, gene_to_LFC, logger)
    # mouse
    gene_to_up_or_down, gene_to_LFC = parse_file(in2, gene_to_up_or_down, gene_to_LFC, logger)
    # bat                                           
    gene_to_up_or_down, gene_to_LFC = parse_file(in3, gene_to_up_or_down, gene_to_LFC, logger)  
    
    parse_orthofinder(orthofinder, outfile, hu_gene_to_prot,
                      mo_gene_to_prot,
                      gene_to_up_or_down, gene_to_LFC,
                      logger)
            
