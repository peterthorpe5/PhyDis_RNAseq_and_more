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


def parse_cluster(in1, in2, in3):
    """fucn take in the  GOI_set
    if the gene in the orthofinder output check if gene in cluster. """
    f_in = open(orthof, "r")
    f_out = open(outfile, "w")
    matrix_found_count = 0
    for line in f_in:
        if test_line(line):
            genes = line.split()
            DE_gene = genes[0]
            OG = genes[1]
            print(OG)
            
                    



if "-v" in sys.argv or "--version" in sys.argv:
    print("v0.0.1")
    sys.exit(0)


usage = """Use as follows:

$ python .py -h

"""

parser = OptionParser(usage=usage)

parser.add_option("--in1", dest="in1",
                  default="Orthogroups.txt",
                  help="Orthogroups.txt")

parser.add_option("--in2", dest="in2",
                  default="Orthogroups.txt",
                  help="Orthogroups.txt")

parser.add_option("--in3", dest="in3",
                  default="Orthogroups.txt",
                  help="Orthogroups.txt")

                  
parser.add_option("-o", dest="outfile",
                  default=None,
                  help="cluster_membership_summary")



(options, args) = parser.parse_args()

in1 = options.in1
in2 = options.in2
in3 = options.in3

logfile = "cluster_membership_counts.log" 
# Run as script
if __name__ == '__main__':
    start_time = time.time()
    # Set up logging
    logger = logging.getLogger('cluster_membership.py: %s'
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
    
    parse_cluster(in1, in2, in3, logger)
            
