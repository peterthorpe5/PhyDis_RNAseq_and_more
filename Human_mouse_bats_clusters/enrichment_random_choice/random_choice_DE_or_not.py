# script to parse and filter DE results for sig FC and FDR
# author Pete Thorpe 2021 Oct
# imports

import os
import sys
from optparse import OptionParser
import datetime
import logging
import logging.handlers
import time
import pandas as pd
import random
from collections import defaultdict
import numpy

def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    return line





if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python filter....py --fdr FDR_threshold --log logFC threshold

"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="rnaseq_results",
                  default="phydis_LOGFC_1.50_FDR_0.01",
                  help="rnaseq_results.txt")

parser.add_option("-l", "--logfc", dest="out",
                  default=1.5,
                  help="logFC threshold: default => 1.5",
                  metavar="FILE")



(options, args) = parser.parse_args()

LOGFC = options.out
rnaseq_results = options.rnaseq_results


logfile = "tmp.log" 
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('random_choice.py: %s'
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
    f_in = open("Genes227.txt", "r")
    GOI_in = []
    for line in f_in:
        line = line.upper()
        GOI_in.append(line.strip())

    f_in = open("all_gene.txt", "r")
    all_genes = []
    for line in f_in:
        line = line.upper()
        all_genes.append(line.strip())

    random_set = set([])
    random_dict = defaultdict(set)
    for j in range(0,100):
        for i in range(0,len(GOI_in)):
            gene = random.choice(all_genes)
            random_dict[j].add(gene.rstrip())

    rnaseq = open(rnaseq_results, "r")
    DE_set = set([])
    for line in rnaseq:
        data = line.split()
        gene = data[0]
        DE_set.add(gene.upper())
    #print(DE_set)
    # interate through th dict
    random_gene_DE = defaultdict(list)
    random_gene_DE_counter = defaultdict(int)
    for random, genes in random_dict.items():
        #print(random, genes)
        for gene in genes:
            #print(gene)
            if gene in DE_set:
                random_gene_DE[random].append(1)
                random_gene_DE_counter[random] += 1
            else:
                random_gene_DE[random].append(0)
    random_de_count = []
    for key, vals in random_gene_DE_counter.items():
        random_de_count.append(vals)
    #print(random_de_count)
    standard_dev = numpy.std(random_de_count)
    the_mean = numpy.mean(random_de_count)
    print("%d random genes obtained 100 times: Of these genes were DE on average (mean): " % len(GOI_in))
    print("standard dev  = ", standard_dev)
    print("mean = ", the_mean)

    count = 0
    for gene in GOI_in:
        if gene.rstrip() in DE_set:
            count = count + 1

    print("from the original 198 GOI list, %d were DE" % count)
    
                
            
    


            
