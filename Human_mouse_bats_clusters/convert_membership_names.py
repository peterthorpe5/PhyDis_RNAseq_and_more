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


def parse_infile(in1, hu_prot_id_to_gene, hu_gene_to_product, 
                 hu_prot_id_to_product, mo_prot_id_to_gene,
                 mo_gene_to_product, 
                 mo_prot_id_to_product,
                 outfile, GOI_literature, logger):
    """ funk to parse the infile and convert the prot id
    to gene names. and annotations. """
    f_in = open(in1, "r")
    f_out = open(outfile, "w")
    for line in f_in:
        if test_line(line):
            if not line.startswith("human\t") and not line.startswith("mouse\t"):
                f_out.write(line)
                #print(line)
                
            else:
                # here we could have any number of species. 
                new_line = line
                data = line.split("\t")
                if len(data) > 2:
                    if data[0] == "human":
                        hu_gene = hu_prot_id_to_gene[data[1].rstrip()]
                        GOI_found = "\n\tGENE_OF_INTEREST_FROM_LIT\t"
                        if hu_gene.upper() in GOI_literature:
                            GOI_found = GOI_found + hu_gene
                        if GOI_found == "\n\tGENE_OF_INTEREST_FROM_LIT\t":
                            GOI_found =  ""
                        hu_annot = hu_prot_id_to_product[data[1].rstrip()]
                        # print(data[1], hu_gene, hu_annot)
                        outfmt = "\t".join([data[1].rstrip(), hu_gene, hu_annot])
                        new_line = new_line.replace(data[1].rstrip(), outfmt)
                        new_line =  new_line.rstrip()  + GOI_found + "\n"

                    if data[2] == "mouse":
                        mo_gene = mo_prot_id_to_gene[data[3].rstrip()]
                        GOI_found = "\n\tGENE_OF_INTEREST_FROM_LIT\t"
                        if mo_gene.upper() in GOI_literature:
                            GOI_found = GOI_found + mo_gene
                        if GOI_found == "\n\tGENE_OF_INTEREST_FROM_LIT\t":
                            GOI_found =  ""
                        mo_annot = mo_prot_id_to_product[data[3].rstrip()]
                        # print(data[3], mo_gene, mo_annot)
                        outfmt = "\t".join([data[3].rstrip(), mo_gene, mo_annot])
                        new_line = new_line.replace(data[3].rstrip(), outfmt)
                        new_line =  new_line.rstrip()  + GOI_found + "\n"
                    
                    if data[0] == "mouse":
                        mo_gene = mo_prot_id_to_gene[data[1].rstrip()]
                        GOI_found = "\n\tGENE_OF_INTEREST_FROM_LIT\t"
                        if mo_gene.upper() in GOI_literature:
                            GOI_found = GOI_found + mo_gene
                        if GOI_found == "\n\tGENE_OF_INTEREST_FROM_LIT\t":
                            GOI_found =  ""
                        mo_annot = mo_prot_id_to_product[data[1].rstrip()]
                        # print(data[1], mo_gene, mo_annot)
                        outfmt = "\t".join([data[1].rstrip(), mo_gene, mo_annot])
                        new_line = new_line.replace(data[1].rstrip(), outfmt)
                        new_line =  new_line.rstrip()  + GOI_found + "\n"
                    # print(new_line.rstrip())
                    f_out.write(new_line)
                else:
                    # print(line)
                    if data[0] == "human":
                        hu_gene = hu_prot_id_to_gene[data[1].rstrip()]
                        GOI_found = "\n\tGENE_OF_INTEREST_FROM_LIT\t"
                        if hu_gene.upper() in GOI_literature:
                            GOI_found = GOI_found + hu_gene
                        if GOI_found == "\n\tGENE_OF_INTEREST_FROM_LIT\t":
                            GOI_found =  ""
                        hu_annot = hu_prot_id_to_product[data[1].rstrip()]
                        # print(data[1], hu_gene, hu_annot)
                        outfmt = "\t".join([data[1].rstrip(), hu_gene, hu_annot])
                        new_line = new_line.replace(data[1].rstrip(), outfmt)
                        new_line =  new_line.rstrip()  + GOI_found + "\n"

                    if data[0] == "mouse":
                        mo_gene = mo_prot_id_to_gene[data[1].rstrip()]
                        GOI_found = "\n\tGENE_OF_INTEREST_FROM_LIT\t"
                        if mo_gene.upper() in GOI_literature:
                            GOI_found = GOI_found + mo_gene
                        if GOI_found == "\n\tGENE_OF_INTEREST_FROM_LIT\t":
                            GOI_found =  ""
                        mo_annot = mo_prot_id_to_product[data[1].rstrip()]
                        # print(data[1], mo_gene, mo_annot)
                        outfmt = "\t".join([data[1].rstrip(), mo_gene, mo_annot])
                        new_line = new_line.replace(data[1].rstrip(), outfmt)
                        new_line =  new_line.rstrip()  + GOI_found + "\n"

                    # print(new_line.rstrip())
                    new_line = new_line.replace("\n\n", "\n")
                    f_out.write(new_line)

                    
        
            
        
    
    

if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python .py -h

"""

parser = OptionParser(usage=usage)

parser.add_option("--in", dest="in1",
                  default=None,
                  help="Orthogroups.txt")

parser.add_option("--hu", dest="human",
                  default="human_gene_to_symblo.info",
                  help="human_gene_to_symblo.info")

parser.add_option("--mu", dest="mouse",
                  default="mus_gene_to_symblo.info",
                  help="mus_gene_to_symblo.info")

                  
parser.add_option("-o", dest="outfile",
                  default="test",
                  help="cluster_membership_summary")



(options, args) = parser.parse_args()

in1 = options.in1
human = options.human
mouse = options.mouse
outfile = options.outfile

logfile = "convert_membership_names.log" 
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('convert_membership_names.py: %s'
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
    GOI_literature = set([])
    f_goi = open("Genes227.txt", "r")
    for line in f_goi:
        if test_line(line):
            line = line.upper()
            GOI_literature.add(line.strip())
    
    hu_gene_to_prot, hu_prot_id_to_gene, hu_gene_to_product, \
            hu_prot_id_to_product = parse_NCBI_gffcol9_info(human)

    mo_gene_to_prot, mo_prot_id_to_gene, mo_gene_to_product, \
            mo_prot_id_to_product = parse_NCBI_gffcol9_info(mouse)
            
    parse_infile(in1, hu_prot_id_to_gene, hu_gene_to_product, 
                 hu_prot_id_to_product, mo_prot_id_to_gene,
                 mo_gene_to_product, 
                 mo_prot_id_to_product,
                 outfile, GOI_literature, logger)
