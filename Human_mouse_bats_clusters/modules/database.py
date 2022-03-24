#!/usr/bin/env python3
#
# metapy_tools.py
#
#
# 2022
# Author:  Peter Thorpe

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



def parse_NCBI_gffcol9_info(infile):
    """func to get the required info from a preparpared info file, as such

    zcat GCF_000001635.27_GRCm39_genomic.gff.gz | grep "CDS"
    | cut -f9 | uniq > mus_gene_to_symblo.info
    
    This will result sin output like this:
    ID=cds-XP_006503333.1;Parent=rna-XM_006503270.5;Dbxref=GeneID:66260,
    Genbank:XP_006503333.1,MGI:MGI:1913510;Name=XP_006503333.1;gbkey=CDS;
    gene=Tmem54;product=transmembrane protein 54 isoform X1;
    protein_id=XP_006503333.1

    we are interested in saving the protein_id= cds- gene= Name= product=

    The AA file does not match the gene ids used by the group, so
    this info is required for a traslation.
    e.g. human "A3GALT2" results in this:
    
    ID=cds-NP_001073907.1;Parent=rna-NM_001080438.1;Dbxref=CCDS:CCDS60080.1,
    Ensembl:ENSP00000475261.1,GeneID:127550,Genbank:NP_001073907.1,
    HGNC:HGNC:30005;Name=NP_001073907.1;gbkey=CDS;gene=A3GALT2;
    product=alpha-1%2C3-galactosyltransferase 2;protein_id=NP_001073907.1;
    tag=MANE Select
    """
    # set up some default dicts
    gene_to_prot = defaultdict(str)
    prot_id_to_gene = defaultdict(int)
    gene_to_product = defaultdict(str)
    prot_id_to_product = defaultdict(int)
    
    f_in = open(infile, "r")
    input_count = 0
    goi_set = set([])
    for line in f_in:
        input_count = input_count + 1
        if test_line(line):
            try:
                cds_info = line.split("ID=cds-")[1]
                cds_info = cds_info.split(";")[0]
            except:
                # no cds info
                continue
            try:
                # no gene info
                gene = line.split("gene=")[1]
                gene = gene.split(";")[0]
            except:
                continue
            try:
                product = line.split("product=")[1]
                product = product.split(";")[0]
            except:
                # no product info
                product = ""
            # print("cds_info, gene, product")
            # print(cds_info, gene, product)
            gene_to_prot[gene] = cds_info
            prot_id_to_gene[cds_info] = gene
            gene_to_product[gene] = product
            prot_id_to_product[cds_info] = product
    return gene_to_prot, prot_id_to_gene, gene_to_product, prot_id_to_product

