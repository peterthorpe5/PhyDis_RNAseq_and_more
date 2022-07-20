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


def test_line(line):
    """returns true lines. Not comments or blank line"""
    if not line.strip():
        return False  # if the last line is blank
    if line.startswith("#"):
        return False  # comment line
    return line


def parse_file(infile, outfile, species, min_tpm=5):
    """func to parse the de results and writes the sig results"""
    out = open(outfile, "w")
    df = pd.read_table(infile)
    df = df.reset_index()
    f_out = open(outfile, "w")
    gene_names = set([])
    species_cortex = species + "_ctx"
    species_striatum = species + "_str"
    count = 0
    for index, row in df.iterrows():
        count = count + 1
        #print(index, row)
        # incositant naming in files, so use pandas
        gene_name = row["Name"]
        sp_cortex = row["%s_ctx" % species]
        sp_striatum = row["%s_str" % species]
        if int(sp_cortex) + int(sp_striatum) > int(min_tpm):
            out_data = "%s\t%s\t%s\n" % (gene_name, sp_cortex, sp_striatum)
            #print(out_data)
            gene_names.add(gene_name)
    for gene in gene_names:
        f_out.write( gene + "\n")
    f_out.close()
    return gene_names, count
        
    

human, count = parse_file("tpm_matrix_14k.txt", "human_head_expressed", "Human", 5)
Mouse, count = parse_file("tpm_matrix_14k.txt", "Mouse_head_expressed", "Mouse", 5)
PhyDis, count = parse_file("tpm_matrix_14k.txt", "PhyDis_head_expressed", "PhyDis", 5)

print( "out of %d staring gene list" % count)
print("%d genes expressed in mouse head" % len(Mouse))
print("%d genes expressed in human head" % len(human))
print("%d genes expressed in PhyDis head" % len(PhyDis))



