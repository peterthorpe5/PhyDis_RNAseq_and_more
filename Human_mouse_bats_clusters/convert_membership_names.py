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
    
    hu_gene_to_prot, hu_prot_id_to_gene, hu_gene_to_product, \
            hu_prot_id_to_product = parse_NCBI_gffcol9_info(human)

    mo_gene_to_prot, mo_prot_id_to_gene, mo_gene_to_product, \
            mo_prot_id_to_product = parse_NCBI_gffcol9_info(mouse)   
