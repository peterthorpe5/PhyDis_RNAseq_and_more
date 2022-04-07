set -e
# set up some dir and variables if required:
# Notes for Paolo. Put the raw RNAseq output into a working directory.
# Then run the command as below. The FDR and LOGFC can be set up in a loop
# so a whole range of thresholds can be tested. 
# What do we consider DE??

values="0.40 0.60 0.80 1.00 1.20 1.40 1.60 1.80 2.00"
FDR=0.01
RNAseq_dir=/storage/home/users/sonia_vernes/human_mouse_bat_clustering/filtered_by_logFC_${LOGFC}_FDR_${FDR}/

# pandas is required
#conda activate python36

for LOGFC in ${values}
	do
    RNAseq_dir=/storage/home/users/sonia_vernes/human_mouse_bat_clustering/filtered_by_logFC_${LOGFC}_FDR_${FDR}/
    echo "RNAseq_dir"
    echo "python filter_RNAseq_CSV.py --fdr ${FDR} -l ${LOGFC}"
    python filter_RNAseq_CSV.py --fdr ${FDR} -l ${LOGFC}

    # 4)
    # convert the "gene" name to protien Id.  - note some names dont convert.
    # the wonders of data!

    echo "python convert_gene_names.py -i human_gene_to_symblo.info \
        --goi ${RNAseq_dir}/human_LOGFC_${LOGFC}_FDR_${FDR} \
        -o ${RNAseq_dir}/human_prot_id_LOGFC_${LOGFC}_FDR_${FDR}"
    python convert_gene_names.py -i human_gene_to_symblo.info \
        --goi ${RNAseq_dir}/human_LOGFC_${LOGFC}_FDR_${FDR} \
        -o ${RNAseq_dir}/human_prot_id_LOGFC_${LOGFC}_FDR_${FDR}
        
    # repeat for mouse
    echo "python convert_gene_names.py -i mus_gene_to_symblo.info   \  
        --goi ${RNAseq_dir}/mouse_LOGFC_${LOGFC}_FDR_${FDR} \
        -o ${RNAseq_dir}/mouse_prot_id_LOGFC_${LOGFC}_FDR_${FDR}"
    python convert_gene_names.py -i mus_gene_to_symblo.info   --goi ${RNAseq_dir}/mouse_LOGFC_${LOGFC}_FDR_${FDR} \
        -o ${RNAseq_dir}/mouse_prot_id_LOGFC_${LOGFC}_FDR_${FDR}

# repeat for bats. 
# no need for bats. Already in correct format. 

# 5)
# Now we have prot ID names that have been filtered on thresholds, which are viable. Lets pick the clusters 
# which have DE genes in them. 
    python interogate_clusters.py --goi ${RNAseq_dir}/human_prot_id_LOGFC_${LOGFC}_FDR_${FDR} \
        -o human_DE_${LOGFC}_FDR_${FDR}_clusters

    python interogate_clusters.py --goi ${RNAseq_dir}/mouse_prot_id_LOGFC_${LOGFC}_FDR_${FDR} \
        -o mouse_DE_${LOGFC}_FDR_${FDR}_clusters

    python interogate_clusters.py --goi ${RNAseq_dir}/phydis_LOGFC_${LOGFC}_FDR_${FDR} \
        -o phydis_DE_${LOGFC}_FDR_${FDR}_clusters

# 6)
# Now we have to clusters, lets see which species could be found in these clusters. 

    python cluster_membership.py --in1 human_DE_${LOGFC}_FDR_${FDR}_clusters \
        --in2 mouse_DE_${LOGFC}_FDR_${FDR}_clusters \
        --in3 phydis_DE_${LOGFC}_FDR_${FDR}_clusters \
        -o cluster_summary_${LOGFC}_FDR_${FDR}.RESULTS
        
    # 7) optional: convert the protein IDS for humans and mouse back to "gene name"

    python convert_membership_names.py --hu human_gene_to_symblo.info \
        --mu mus_gene_to_symblo.info \
        --in cluster_summary_${LOGFC}_FDR_${FDR}.RESULTS \
        -o cluster_summary_${LOGFC}_FDR_${FDR}_gene_name.RESULTS
   
   # 8) add the direction of DE
   python convert_membership_names_up_or_down.py  --infile cluster_summary_${LOGFC}_FDR_${FDR}_gene_name.RESULTS  \
  --in1 ${RNAseq_dir}/human_prot_id_LOGFC_${LOGFC}_FDR_${FDR} \
    --in2  ${RNAseq_dir}/mouse_prot_id_LOGFC_${LOGFC}_FDR_${FDR}  \
    --in3 ${RNAseq_dir}/phydis_LOGFC_${LOGFC}_FDR_${FDR} \
    -o cluster_summary_${LOGFC}_FDR_${FDR}_gene_name_up_or_down.RESULTS
done



