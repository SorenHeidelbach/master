#!/bin/bash
wd="/shared-nfs/SH/nanodisco/MGM1"
cd $wd

bash /shared-nfs/SH/code/01_bin_quality_check.sh -i $wd/t_SNE_clusters_10kb_15nb -o $wd/t_SNE_clusters_10kb_15nb -x fasta
