#! -cwd
cd /storage/home/users/sonia_vernes/human_mouse_bat_clustering/

conda activate orthofinder

orthofinder -S diamond -t 24 -f /storage/home/users/sonia_vernes/human_mouse_bat_clustering/AA/

 python ../mcl_to_cafe.py -sp "GROS Hsc GPLIN GPALN Minc Hetgly" -i ./Results_*/Orthogroups.txt -o nematodes_orthofinder.clusters.txt
 