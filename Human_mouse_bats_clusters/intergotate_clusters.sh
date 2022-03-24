# set up some dir and variables if required:
# Notes for Paolo. Put the raw RNAseq output into a working directory.
# Then run the command as below. The FDR and LOGFC can be set up in a loop
# so a whole range of thresholds can be tested. 
# What do we consider DE??


FDR=0.01
LOGFC=1.50
RNAseq_dir=/storage/home/users/sonia_vernes/human_mouse_bat_clustering/filtered_by_logFC_${LOGFC}_FDR_${FDR}/

# pandas is required
conda activate python36

# 1)
# first run orthofinder on the human, mouse GCF_000001635.27_GRCm39_genomic.gff.gz
# GCF_000001405.39_GRCh38.p13_genomic.gff.gz and Bat PhyDis datasets. 
# These are the versions used. 

# 2)
# get the prot id to gene name infor from the gff for mouse and human:
# zcat GCF_000001635.27_GRCm39_genomic.gff.gz | grep "CDS" | cut -f9 | uniq > mus_gene_to_symblo.info

# 3)
# filter the RNAseq results. The exact thresholds are a work on progress right now. 
# SO scripted so this can be altered:
# RNAseq outfiles end in _allresults.csv
# python filter_RNAseq_CSV.py  (defaults can be changed at the command line)

python filter_RNAseq_CSV.py --fdr ${FDR} -l ${LOGFC}

# 4)
# convert the "gene" name to protien Id.  - note some names dont convert.
# the wonders of data!


python convert_gene_names.py -i human_gene_to_symblo.info --goi ${RNAseq_dir}/human_LOGFC_${LOGFC}_FDR_${FDR} \
    -o ${RNAseq_dir}/human_prot_id_LOGFC_${LOGFC}_FDR_${FDR}
    
# repeat for mouse
python convert_gene_names.py -i mus_gene_to_symblo.info --goi ${RNAseq_dir}/mouse_LOGFC_${LOGFC}_FDR_${FDR} \
    -o ${RNAseq_dir}/mouse_prot_id_LOGFC_${LOGFC}_FDR_${FDR}

# repeat for bats. 
# no need for bats. Already in correct format. 

# 5)
# Now we have prot ID names that have been filtered on thresholds, which are viable. Lets pick the clusters 
# which have DE genes in them. 
python interogate_clusters.py --goi ${RNAseq_dir}/human_prot_id_LOGFC_${LOGFC}_FDR_${FDR} -o human_DE_${LOGFC}_FDR_${FDR}_clusters

python interogate_clusters.py --goi ${RNAseq_dir}/mouse_prot_id_LOGFC_${LOGFC}_FDR_${FDR} -o mouse_DE_${LOGFC}_FDR_${FDR}_clusters

python interogate_clusters.py --goi ${RNAseq_dir}/phydis_LOGFC_${LOGFC}_FDR_${FDR} -o phydis_DE_${LOGFC}_FDR_${FDR}_clusters

# 6)
# Now we have to clusters, lets see which species could be found in these clusters. 

python cluster_membership.py --in1 human_DE_${LOGFC}_FDR_${FDR}_clusters \
    --in2 mouse_DE_${LOGFC}_FDR_${FDR}_clusters \
    --in3 phydis_DE_${LOGFC}_FDR_${FDR}_clusters



