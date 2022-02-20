



# Paths to fastq annotated fast5 and reference
fast5_nat=/shared-nfs/SH/nanodisco/dataset/EC_NAT/
fast5_wga=/shared-nfs/SH/nanodisco/dataset/EC_WGA/
ref=/shared-nfs/SH/nanodisco/reference/Ecoli_K12_MG1655_ATCC47076.fasta

# Preprocessing
nanodisco preprocess -p 60 -f $fast5_nat -o /shared-nfs/SH/nanodisco/ecoli -r $ref -s ecoli_nat
nanodisco preprocess -p 60 -f $fast5_wga -o /shared-nfs/SH/nanodisco/ecoli -r $ref -s ecoli_wga

# See number of reference chunks
nanodisco chunk_info -r $ref

# Compute current difference between NAT and WGA
nanodisco difference -nj 10 -nc 1 -p 2 -f 1 -l 932 \
  -i /shared-nfs/SH/nanodisco/ecoli \
  -o /shared-nfs/SH/nanodisco/difference \
  -w ecoli_wga \
  -n ecoli_nat \
  -r $ref
  

# Merge current difference  
nanodisco merge -d /shared-nfs/SH/nanodisco/difference -o /shared-nfs/SH/nanodisco/difference_merged -b EC_subset


nanodisco motif -p 40 \
  -b EC_subset \
  -d /shared-nfs/SH/nanodisco/difference_merged/EC_subset_difference.RDS \
  -o /shared-nfs/SH/nanodisco/motif \
  -r $ref -a
  
  
nanodisco characterize -p 40 \
  -b EC_subset \
  -d /shared-nfs/SH/nanodisco/difference_merged/EC_subset_difference.RDS \
  -o /shared-nfs/SH/nanodisco/motif_finemapping \
  -m GATC,CCWGG \
  -t nn,rf,knn \
  -r $ref
  
  
  
  
# =================================================================================================
# = 3. methylation PROFILE

# Coverage
nanodisco coverage -b /shared-nfs/SH/nanodisco/ecoli/ecoli_nat.sorted.bam -r $ref -o /shared-nfs/SH/nanodisco/
nanodisco coverage -b /shared-nfs/SH/nanodisco/ecoli/ecoli_wga.sorted.bam -r $ref -o /shared-nfs/SH/nanodisco/

# Methylation profile 4mer
nanodisco profile -p 4 \
  -r $ref \
  -d /shared-nfs/SH/nanodisco/difference_merged/EC_subset_difference.RDS \
  -w /shared-nfs/SH/nanodisco/ecoli_wga.cov \
  -n /shared-nfs/SH/nanodisco/ecoli_nat.cov \
  -b EC_subset \
  -o /shared-nfs/SH/nanodisco/profile_4mer \
  -a 4mer
  
# Methylation profile 6mer
nanodisco profile -p 4 \
  -r $ref \
  -d /shared-nfs/SH/nanodisco/difference_merged/EC_subset_difference.RDS \
  -w /shared-nfs/SH/nanodisco/ecoli_wga.cov \
  -n /shared-nfs/SH/nanodisco/ecoli_nat.cov \
  -b EC_subset \
  -o /shared-nfs/SH/nanodisco/profile_6mer \
  -a 6mer

# methylation profile all
nanodisco profile -p 50 \
  -r $ref \
  -d /shared-nfs/SH/nanodisco/difference_merged/EC_subset_difference.RDS \
  -w /shared-nfs/SH/nanodisco/ecoli_wga.cov \
  -n /shared-nfs/SH/nanodisco/ecoli_nat.cov \
  -b EC_subset \
  -o /shared-nfs/SH/nanodisco/profile_all \
  -a all


# Scoring each motif loci based on predefined motif
nanodisco score -p 10 -b EC_subset \
  -d /shared-nfs/SH/nanodisco/difference_merged/EC_subset_difference.RDS \
  -o /shared-nfs/SH/nanodisco/scores_GATC \
  -m GATC \
  -r $ref
nanodisco score -p 10 -b EC_subset \
  -d /shared-nfs/SH/nanodisco/difference_merged/EC_subset_difference.RDS \
  -o /shared-nfs/SH/nanodisco/scores_GATC_CCWGG \
  -m GATC,CCWGG \
  -r $ref



# =================================================================================================
# = MGM1 metabat
module load MetaBAT/2.12.1-foss-2018a
metabat -i /shared-nfs/SH/nanodisco/MGM1/reference/metagenome.fasta \
  -o /shared-nfs/SH/nanodisco/MGM1/metabat/metabat \
  -t 20 --saveCls

module load CheckM/1.1.2-foss-2018a-Python-3.6.4
checkm lineage_wf /shared-nfs/SH/nanodisco/MGM1/metabat \
  /shared-nfs/SH/nanodisco/MGM1/checkm \
  -x fa > /shared-nfs/SH/nanodisco/MGM1/checkm/summary.txt















