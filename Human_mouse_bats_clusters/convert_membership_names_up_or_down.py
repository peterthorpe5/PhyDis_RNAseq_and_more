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
import pandas as pd


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


def parse_infile(infile, outfile, gene_to_up_or_down, gene_to_LFC,
                 logger):
    """ funk to parse the infile and convert the prot id
    to gene names. and annotations. """
    f_in = open(infile, "r")
    f_out = open(outfile, "w")
    for line in f_in:
        if test_line(line):
            if not line.startswith("human\t") and not line.startswith("mouse\t") and not line.startswith("phydis\t"):
                f_out.write(line)
                #print(line)
                
            else:
                # here we could have any number of species. 
                new_line = line
                line = line.rstrip()
                data = line.split("\t")
                if data[1] in gene_to_up_or_down:
                    element = data[1].strip()
                    #print(element)
                    if gene_to_up_or_down[element]:
                        #print("yes")
                        direction = gene_to_up_or_down[element]
                        LFC = gene_to_LFC[element]
                        

                        outfmt = "\t".join([element, direction, "LogFC = %.3f" % (LFC)])
                        new_line = new_line.replace(element, outfmt)

                        #logger.info(new_line)
                        #f_out.write(new_line)
                print(data)
                if len(data) < 5: 
                    logger.info(new_line)
                    f_out.write(new_line)
                    continue
                if data[5] in gene_to_up_or_down:
                   element = data[5].strip()
                    #print(element)
                   if gene_to_up_or_down[element]:
                        #print("yes")
                       direction = gene_to_up_or_down[element]
                       LFC = gene_to_LFC[element]
                        

                       outfmt = "\t".join([element, direction, "LogFC = %.3f" % (LFC)])
                       new_line = new_line.replace(element, outfmt)

                       #logger.info(new_line)
                       #f_out.write(new_line)
                
                if len(data) < 9: 
                    logger.info(new_line)
                    f_out.write(new_line)
                    continue
                    
                if data[9] in gene_to_up_or_down:
                    element = data[9].strip()
                    #print(element)
                    if gene_to_up_or_down[element]:
                        #print("yes")
                        direction = gene_to_up_or_down[element]
                        LFC = gene_to_LFC[element]
                        

                        outfmt = "\t".join([element, direction, "LogFC = %.3f" % (LFC)])
                        new_line = new_line.replace(element, outfmt)

                        logger.info(new_line)
                        f_out.write(new_line)
 
               


if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python .py -h

"""

parser = OptionParser(usage=usage)

parser.add_option("--infile", dest="infile",
                  default=None,
                  help="Orthogroups.txt")

parser.add_option("--in1", dest="in1",
                  default=None,
                  help="human filtered RNAseq results")

parser.add_option("--in2", dest="in2",
                  default=None,
                  help="mouse filtered RNAseq results")

parser.add_option("--in3", dest="in3",
                  default=None,
                  help="phydis filtered RNAseq results")

                  
parser.add_option("-o", dest="outfile",
                  default="test",
                  help="cluster_membership_summary")



(options, args) = parser.parse_args()

infile = options.infile
in1 = options.in1
in2 = options.in2
in3 = options.in3
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
    
    
    # set up some dictionaries to capture the data. 
    gene_to_up_or_down = defaultdict(str)
    gene_to_LFC = defaultdict(float)
    # add the human data
    gene_to_up_or_down, gene_to_LFC = parse_file(in1, gene_to_up_or_down, gene_to_LFC, logger)
    # mouse
    gene_to_up_or_down, gene_to_LFC = parse_file(in2, gene_to_up_or_down, gene_to_LFC, logger)
    # bat                                           
    gene_to_up_or_down, gene_to_LFC = parse_file(in3, gene_to_up_or_down, gene_to_LFC, logger)                                                  
                                                 
    parse_infile(infile, outfile, gene_to_up_or_down, gene_to_LFC,
                 logger)
