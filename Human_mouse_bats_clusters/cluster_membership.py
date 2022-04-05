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


def populate_set(infile, set, indict):
    """func to return a set of clusters"""
    f_in = open(infile, "r")
    for line in f_in:
        if test_line(line):
            genes = line.split()
            DE_gene = genes[0]
            OG = genes[1].split(":")[0]
            # print(OG)
            set.add(OG)
            indict[OG] = DE_gene
    return set, indict
            


def parse_DE_cluster(in1, in2, in3, outfile, logger):
    """fucn take in the  GOI_set
    if the gene in the orthofinder output check if gene in cluster. """
    f_out = open(outfile, "w")
    species1 = in1.split("_")[0]
    species2 = in2.split("_")[0]
    species3 = in3.split("_")[0]
    logger.info("species1 = %s " % species1)
    logger.info("species2 = %s " % species2)
    logger.info("species3 = %s " % species3)
    f_out.write("species1 = %s \n" % species1)
    f_out.write("species2 = %s \n" % species2)
    f_out.write("species3 = %s \n" % species3)
    # set up some sets
    species1_set = set([])
    species2_set = set([])
    species3_set = set([])
    # set up some default dicts
    species1_dict = defaultdict()
    species2_dict = defaultdict()
    species3_dict = defaultdict()
    # call the fun to populate the sets
    species1_set, species1_dict = populate_set(in1, species1_set, 
                                               species1_dict)
    species2_set, species2_dict = populate_set(in2, species2_set, 
                                               species2_dict)
    species3_set, species3_dict = populate_set(in3, species3_set, 
                                               species3_dict)
    common_to_all = species1_set.intersection(species2_set, 
                                              species3_set)
    
    outdata = "Clusters which contained a DE genes from ALL species: %d\t%s\n" % (len(common_to_all), common_to_all)                                      
    logger.info(outdata)
    f_out.write(outdata)
    logger.info("these genes are")
    f_out.write("DE genes from ALL species - these genes are\n")
    
    for cluster in common_to_all:
        gene1 = species1_dict[cluster]
        gene2 = species2_dict[cluster]
        gene3 = species3_dict[cluster]
        outdata = "\t".join([species1, gene1, species2, 
                             gene2, species3, gene3])
        logger.info(outdata)
        f_out.write(outdata + "\n")
    
    ###################
    common_to_1_and_2 = species1_set.intersection(species2_set)
    common_to_1_and_2_unique = common_to_1_and_2.difference(common_to_all)
    
    outdata = "Clusters which contained a DE genes from species1 and species2 (BUT NOT THOSE ALREADY IDENTIFIED AS COMMON TO ALL: %d\t%s\n" % (len(common_to_1_and_2_unique), common_to_1_and_2_unique)                                      
    logger.info(outdata)
    f_out.write(outdata)

    logger.info("these genes are")
    f_out.write("\n%s and %s COMMON: these genes are\n" % (species1, species2))
    
    for cluster in common_to_1_and_2_unique:
        gene1 = species1_dict[cluster]
        gene2 = species2_dict[cluster]
        outdata = "\t".join([species1, gene1, species2, 
                             gene2])
        logger.info(outdata)
        f_out.write(outdata + "\n")
    
    ####################
    common_to_1_and_3 = species1_set.intersection(species3_set)
    common_to_1_and_3_unique = common_to_1_and_3.difference(common_to_all)
    
    outdata = "Clusters which contained a DE genes from species1 and species3 (BUT NOT THOSE ALREADY IDENTIFIED AS COMMON TO ALL: %d\t%s\n" % (len(common_to_1_and_3_unique), common_to_1_and_3_unique)                                      
    logger.info(outdata)
    f_out.write(outdata)
    
    logger.info("these genes are")
    f_out.write("\n%s and %s COMMON: these genes are\n" % (species1, species3))
    
    for cluster in common_to_1_and_3_unique:
        gene1 = species1_dict[cluster]
        gene3 = species3_dict[cluster]
        outdata = "\t".join([species1, gene1, species3, 
                             gene3])
        logger.info(outdata)
        f_out.write(outdata + "\n")
    
    ####################
    common_to_2_and_3 = species2_set.intersection(species3_set)
    common_to_2_and_3_unique = common_to_2_and_3.difference(common_to_all)
    
    outdata = "Clusters which contained a DE genes from species2 and species3 (BUT NOT THOSE ALREADY IDENTIFIED AS COMMON TO ALL: %d\t%s\n" % (len(common_to_2_and_3_unique), common_to_2_and_3_unique)                                      
    logger.info(outdata)
    f_out.write(outdata)
    
    logger.info("these genes are")
    f_out.write("\n%s and %s COMMON: these genes are\n" % (species2, species3))
    
    for cluster in common_to_2_and_3_unique:
        gene2 = species2_dict[cluster]
        gene3 = species3_dict[cluster]
        outdata = "\t".join([species2, gene2, species3, 
                             gene3])
        logger.info(outdata)
        f_out.write(outdata + "\n")
    f_out.close()
    

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
                  default="test",
                  help="cluster_membership_summary")



(options, args) = parser.parse_args()

in1 = options.in1
in2 = options.in2
in3 = options.in3
outfile = options.outfile

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
    
    parse_DE_cluster(in1, in2, in3, outfile, logger)
            
