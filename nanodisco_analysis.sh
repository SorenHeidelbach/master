
export threads=20
export wd=/shared-nfs/SH





######################################################################################################
# # Extracting of fastq from fast5
######################################################################################################

# AD WGA reads
python $wd/code/scripts/fast5_to_fastq.py /shared-nfs/AAU_CS/2021-10-22-SH-WGA/AD/20211022_1448_MN24067_FAQ34806_ae2f68b3/recall \ 
  $wd/samples/AD/WGA/20211022_1448_MN24067_FAQ34806_ae2f68b3.fastq \
  --jobs 20 --batch_size 4000

# AD NAT reads
NAT_fast5=$(find /shared-nfs/AAU_CS/211026_NAT_SH_AD_Zymo/AD1/*/recall* -maxdepth 0)
count=0
for i in $NAT_fast5; do 
  fastq=$wd/samples/AD/NAT/$(basename $(dirname $i))_$count
  echo $fastq
  python $wd/code/scripts/fast5_to_fastq.py $i $fastq --jobs 20 --batch_size 4000
  (( count++ ))
done

# Zymo WGA reads
python $wd/code/scripts/fast5_to_fastq.py /shared-nfs/AAU_CS/2021-10-22-SH-WGA/ZYMO/20211022_1448_MN35936_FAQ30346_fe5e66e2/recall \ 
  $wd/samples/zymo/WGA/20211022_1448_MN35936_FAQ30346_fe5e66e2.fastq \
  --jobs 20 --batch_size 4000

# Zymo NAT reads
NAT_fast5=$(find /shared-nfs/AAU_CS/211026_NAT_SH_AD_Zymo/Zymo/recall/*/fast5 -maxdepth 0)
count=0
for i in $NAT_fast5; do 
  fastq=$wd/samples/zymo/WGA/$(basename $(dirname $i))_$count
  echo $fastq
  python $wd/code/scripts/fast5_to_fastq.py $i $fastq --jobs 20 --batch_size 4000
  (( count++ ))
done




######################################################################################################
# Check Zymo Coverage to Reference
######################################################################################################

#cd $wd/samples/zymo/
#wget http://nanopore.s3.climb.ac.uk/mockcommunity/v2/Zymo-Isolates-SPAdes-Illumina.fasta

get_cov () {
  outdir="$(dirname $2)/cov_$3"
  bam_out="$outdir/aln_sorted.bam"
  cov_out="$outdir/reference_coverage.txt"
  map_read_length_out="$outdir/mapped_read_length.tsv"
  log() {
    echo $1 | tee -a $outdir/log.txt
  }
  
  mkdir -p $outdir
  log "REF  $1"
  log "FASTQ  $2"
  log "BAM  $bam_out"
  log "COV  $cov_out"
  # Map reads
  module load Minimap2/2.17-foss-2020b
  module load SAMtools/1.11-foss-2020b
  minimap2 -ax map-ont -t $threads $1 \
    $2 | \
  samtools view --threads $threads -Sb | samtools sort --threads $threads -o $bam_out
  # Extract read length of mapped reads
  samtools view "$bam_out" | awk '{print $3"\t"length($10)}' >  $map_read_length_out
  # Get coverage of mapped references
  module load BEDTools/2.27.1-GCCcore-6.4.0
  bedtools genomecov -ibam $bam_out -d > $cov_out 
  module purge
}

get_cov "$wd/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta" "$wd/samples/zymo/NAT/pooled_reads.fastq" loman_ref
get_cov "$wd/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta" "$wd/samples/zymo/NAT/pooled_reads_100x.fastq" 100x_loman_ref
get_cov "$wd/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta" "$wd/samples/zymo/NAT/pooled_reads_filt_100x_trim.fastq" filt_100x_trim_loman_ref

get_cov "$wd/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta" "$wd/samples/zymo/WGA/pooled_reads.fastq" loman_ref
get_cov "$wd/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta" "$wd/samples/zymo/WGA/pooled_reads_filt_trim.fastq" filt_trim_loman_ref

get_cov "$wd/samples/zymo/headers.reduced_zymo.fna" "$wd/samples/zymo/NAT/pooled_reads.fastq" real_ref
get_cov "$wd/samples/zymo/headers.reduced_zymo.fna" "$wd/samples/zymo/WGA/pooled_reads.fastq" real_ref

# Investigate read length of zymo NAt and zymo WGA
samtools view pooled_reads_filt_trim.ref.sorted.bam | head | awk '{print $10}'


######################################################################################################
# Subset and filt
######################################################################################################
# Zymo subset to 100x of the NAT reads, genome size = 61.944 Mb, or 36.11067 Mb correcting for the fact that fungi have 1/6 of the coverage
# conda activate py38 # installed in the general enviroment
# Subset the zymo sequence data, as 300x coverage i way too much, subset is done to 100x
# rasusa --coverage 100 --genome-size 36MB \
  # -i "$wd/samples/zymo/NAT/zymo_NAT.fastq" \
  # -o "$wd/samples/zymo/NAT/zymo_NAT_100x.fastq"
# get_cov "$wd/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta" "$wd/samples/zymo/NAT/zymo_NAT_100x.fastq" "100x_loman_ref"

# QC of reads with Nanoplot
fun_nanoplot() {
echo IN $1 - OUT $(dirname $1)/nanoplot_$(basename ${1%.*})
  module load NanoPlot/1.38.0-foss-2020b
  NanoPlot --raw --tsv_stats -t $threads --fastq $1 -o $(dirname $1)/nanoplot_$(basename ${1%.*})
  module purge
}
# fun_nanoplot "$wd/samples/zymo/NAT/pooled_reads_100x.fastq"
fun_nanoplot "$wd/samples/zymo/NAT/pooled_reads.fastq"
fun_nanoplot "$wd/samples/zymo/WGA/pooled_reads.fastq"
fun_nanoplot "$wd/samples/AD/NAT/pooled_reads.fastq"
fun_nanoplot "$wd/samples/AD/WGA/pooled_reads.fastq"

# Filtration of read with Nanofilt
fun_nanofilt () {
  fastq_filt=${1%.*}_filt.fastq
  module load nanofilt/2.6.0-foss-2020b-Python-3.8.6
  NanoFilt -q 7 -l 1000 $1 > $fastq_filt
  module purge
  fun_nanoplot $fastq_filt
}
fun_nanofilt "$wd/samples/zymo/NAT/pooled_reads.fastq"
fun_nanofilt "$wd/samples/zymo/WGA/pooled_reads.fastq"
fun_nanofilt "$wd/samples/AD/NAT/pooled_reads.fastq"
fun_nanofilt "$wd/samples/AD/WGA/pooled_reads.fastq"

# Activate appropiate enviroment
conda activate py38 
# Subset the zymo sequence data, as 300x coverage is way too much, subset is done to 100x
rasusa --coverage 100 --genome-size 36MB \
  -i "$wd/samples/zymo/NAT/pooled_reads_filt.fastq" \
  -o "$wd/samples/zymo/NAT/pooled_reads_filt_100x.fastq"
# Calculate coverage on zymo reference to compare with 100x no filtration subset
get_cov "$wd/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta" "$wd/samples/zymo/NAT/pooled_reads_filt_100x.fastq" "filt_100x_loman_ref"




######################################################################################################
# Adapter trimming
######################################################################################################

# Remove adapters from Nanopore
fun_cutadapt () {
  module load cutadapt/3.4-GCCcore-10.2.0-Python-3.8.6
  echo IN"  "$1
  echo OUT" "${1%.*}_trim.fastq
  cutadapt -g "XGTGAATGTACTTCGTTCAGTTACGTATTGCT;min_overlap=10" -a "GCAATACGTAACTGAACGAAGTX" -e 0.2 -j 20 -o ${1%.*}_trim.fastq $1 > $(dirname $1)/trimming_summary.txt
  module purge
}
fun_cutadapt "$wd/samples/zymo/NAT/pooled_reads_filt_100x.fastq"
get_cov "$wd/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta" "$wd/samples/zymo/NAT/pooled_reads_filt_100x_trim.fastq" "filt_100x_trim_loman_ref"
fun_cutadapt "$wd/samples/zymo/WGA/pooled_reads_filt.fastq"
get_cov "$wd/samples/zymo/Zymo-Isolates-SPAdes-Illumina.fasta" "$wd/samples/zymo/WGA/pooled_reads_filt_100x_trim.fastq" "filt_100x_trim_loman_ref"

fun_cutadapt "$wd/samples/AD/NAT/pooled_reads_filt.fastq"
fun_cutadapt "$wd/samples/AD/WGA/pooled_reads_filt.fastq"




######################################################################################################
# Assembly
######################################################################################################
dos2unix $wd/code/scripts/assembly_and_polish.sh # To make sure script is in right format

## AD WGA
#bash $wd/code/scripts/assembly_and_polish.sh -q "$wd/samples/AD/WGA/pooled_reads_filt.fastq" \
#  -o "$wd/samples/AD/WGA/assembly"
#  -t $threads

# AD NAT
bash $wd/code/scripts/assembly_and_polish.sh -q "$wd/samples/AD/NAT/pooled_reads_filt_trim.fastq" \
  -o "$wd/samples/AD/NAT/assembly_filt_trim"
  -t 30
# Investigate coverage of both WGA and NAT on assembly
get_cov "$wd/samples/AD/NAT/assembly_filt_trim/flye/assembly.fasta" "$wd/samples/AD/NAT/pooled_reads_filt_trim.fastq" "assembly"
get_cov "$wd/samples/AD/NAT/assembly_filt_trim/flye/assembly.fasta" "$wd/samples/AD/WGA/pooled_reads_filt_trim.fastq" "assembly"


# Zymo WGA
# bash $wd/code/scripts/assembly_and_polish.sh -q "$wd/samples/zymo/zymo_WGA_filt.fastq" \
 # -o "$wd/samples/zymo/assembly_WGA"
 # -t $threads

# Zymo NAT
# bash $wd/code/scripts/assembly_and_polish.sh -q "$wd/samples/zymo/pooled_reads_filt.fastq" \
  # -o "$wd/samples/zymo/assembly_NAT"
  # -t $threads

# Zymo NAT 100x
bash $wd/code/scripts/assembly_and_polish.sh -q "$wd/samples/zymo/NAT/pooled_reads_filt_100x_trim.fastq" \
  -o "$wd/samples/zymo/NAT/assembly_filt_100x_trim" \
  -t 80




######################################################################################################
# Prepare fast5s
######################################################################################################

# Get read IDs
extract_ID () {
  echo $1
  awk '/^\@/ {print $1}' $1 | \
  cut -d "_" -f 1 | \
  cut -d "@" -f 2 >  ${1%.*}.ID.txt
}
extract_ID "$wd/samples/zymo/NAT/pooled_reads_filt_100x.fastq"
extract_ID "$wd/samples/zymo/WGA/pooled_reads_filt.fastq"
extract_ID "$wd/samples/AD/NAT/pooled_reads_filt.fastq"
extract_ID "$wd/samples/AD/WGA/pooled_reads_filt.fastq"

# Subset fast5
subset_fast5 () {
  module purge
  module load ont_fast5_api/1.4.7-foss-2018a-Python-3.6.4
  name=$(basename $2)
  echo FROM $1
  echo BY   $2
  echo TO   $(dirname $2)/fast5_${name%.*.*}
  mkdir -p $(dirname $2)/fast5_${name%.*.*}
  fast5_subset -i $1 -s $(dirname $2)/fast5_${name%.*.*} -f $(basename $1)_ -l $2
}
# Export function to use with xargs
export -f subset_fast5

# zymo NAT
find /shared-nfs/AAU_CS/211026_NAT_SH_AD_Zymo/Zymo/recall/*/fast5 -maxdepth 0 | \
  xargs -I {} bash -c "subset_fast5 {} $wd/samples/zymo/NAT/pooled_reads_filt_100x.ID.txt"
# zymo WGA
subset_fast5 /shared-nfs/AAU_CS/2021-10-22-SH-WGA/ZYMO/20211022_1448_MN35936_FAQ30346_fe5e66e2/recall $wd/samples/zymo/WGA/pooled_reads_filt.ID.txt
# AD NAT
find /shared-nfs/AAU_CS/211026_NAT_SH_AD_Zymo/AD1/*/recall* -maxdepth 0 | \
  xargs -I {} bash -c "subset_fast5 {} $wd/samples/AD/NAT/pooled_reads_filt.ID.txt"
# AD WGA
subset_fast5 /shared-nfs/AAU_CS/2021-10-22-SH-WGA/AD/20211022_1448_MN24067_FAQ34806_ae2f68b3/recall $wd/samples/AD/WGA/pooled_reads_filt.ID.txt





######################################################################################################
# Nanopolish Event ALign
######################################################################################################

run_nanopolish(){
  out_bam=${1%.*}.ref.sorted.bam
  out_event_align=$(dirname $3)/$4.txt
  echo "FASTQ $1"
  echo "FAST5 $2"
  echo "REF   $3"
  echo "BAM   $out_bam"
  echo "EVENT $out_event_align"
  echo
  
  # NP Index 
  # echo "  NP index"
  # module load nanopolish/0.11.3-foss-2018a-Python-3.6.4
  # nanopolish index -d $2 $1
  # module purge
  
  # Mapping
  echo "  Mapping to reference"
  module load Minimap2/2.17-foss-2020b
  module load SAMtools/1.11-foss-2020b
  minimap2 -ax map-ont -t $threads $3 $1 | \
    samtools sort -o $out_bam -T $(dirname $3)/read.tmp
  samtools index $out_bam
  
  samtools quickcheck $out_bam
  module purge
  
  # Event align
  echo "  Event align"
  module load nanopolish/0.11.3-foss-2018a-Python-3.6.4
  nanopolish eventalign \
    --reads $1 \
    --bam $out_bam \
    --genome $3 \
    --scale-events > $(dirname $3)/$4.txt
}
# Zymo
## WGA
run_nanopolish "$wd/samples/zymo/WGA/pooled_reads_filt_trim.fastq" \
  "/shared-nfs/AAU_CS/2021-10-22-SH-WGA/ZYMO/20211022_1448_MN35936_FAQ30346_fe5e66e2/fast5" \
  "$wd/samples/zymo/assembly_100x.fasta" \
  "WGA_event_align"
  
time awk '
! ($1 in X) { X[$1] = sprintf ("contig_%.4d.txt", q++ % 100); }
{ print > X[$1]; }' "/shared-nfs/SH/samples/zymo/WGA_event_align.txt"


## NAT
### Copy fast5 to one folder for easier acces (Nanopolish index does not take arrays as input...)
mkdir -p "$wd/samples/zymo/NAT/fast5" 
cp $(echo /shared-nfs/AAU_CS/211026_NAT_SH_AD_Zymo/Zymo/2021*/fast5/*) "$wd/samples/zymo/NAT/fast5"
run_nanopolish "$wd/samples/zymo/NAT/pooled_reads_filt_100x_trim.fastq" \
  "$wd/samples/zymo/NAT/fast5" \
  "$wd/samples/zymo/assembly_100x.fasta" \
  "NAT_event_align"



######################################################################################################
# Nanodisco preprocess
######################################################################################################

# Nanodisco require fast5 to be gzip compressed and not VBZ-compressed
conda activate py38
compress_fast5_fun() {
  compress_fast5 -t 10 --recursive \
    -i $1 \
    -s "$1"_gzip \
    -c gzip
}

compress_fast5_fun "/shared-nfs/SH/samples/zymo/NAT/fast5_multi_100x_filt/" 
compress_fast5_fun "/shared-nfs/SH/samples/zymo/WGA/fast5_multi_filt/"

# Move Fast5 folder to directionary with faster I/O
# mv /shared-nfs/SH/samples/zymo/WGA/fast5_pooled_reads_filt /scratch/users/sheide17/fast5_pooled_reads_filt
# mv /shared-nfs/SH/samples/zymo/NAT/fast5_pooled_reads_filt_100x  /scratch/users/sheide17/zymo_NAT_fast5_pooled_reads_filt_100x 


# Make sure Singularity is in PATH
echo 'export PATH=/shared-nfs/SH/software/singularity/bin:$PATH' >> ~/.bashrc && \
source ~/.bashrc

# Initiate Nanodisco container
singularity run --no-home -B /:/home/nanodisco/dataset -w /shared-nfs/SH/software/nanodisco/nd_env # BIO-server
nd_path=/home/nanodisco/dataset


# Preprocess Zymo 
# WGA
nanodisco preprocess -p 40 \
  -f $nd_path/shared-nfs/SH/samples/zymo/WGA/fast5_multi_filt_gzip \
  -o $nd_path/shared-nfs/SH/samples/zymo/nanodisco/preprocess \
  -r $nd_path/shared-nfs/SH/samples/zymo/nanodisco/assembly_100x.fasta \
  -s zymo_wga

# NAT
nanodisco preprocess -p 40 \
  -f $nd_path/shared-nfs/SH/samples/zymo/NAT/fast5_multi_100x_filt_gzip   \
  -o $nd_path/shared-nfs/SH/samples/zymo/nanodisco/preprocess \
  -r $nd_path/shared-nfs/SH/samples/zymo/nanodisco/assembly_100x.fasta \
  -s zymo_nat

######################################################################################################
# Nanodisco preprocess (on CLAAUDIA AI cloud)
######################################################################################################

# Move NAT fast5, WGA fast5 and assembly
# mkdir -p /user/student.aau.dk/sheide17/projects/modification_binning/zymo/NAT/fast5_basecalled
# scp -r sheide17@student.aau.dk@172.19.20.204:/scratch/users/sheide17/zymo_NAT_fast5_pooled_reads_filt_100x \
  # /user/student.aau.dk/sheide17/projects/modification_binning/zymo/NAT/fast5_basecalled
  
# mkdir -p /user/student.aau.dk/sheide17/projects/modification_binning/zymo/WGA/fast5_basecalled
# scp -r sheide17@student.aau.dk@172.19.20.204:/scratch/users/sheide17/fast5_pooled_reads_filt \
  # /user/student.aau.dk/sheide17/projects/modification_binning/zymo/WGA/fast5_basecalled



# singularity run --no-home -B /user/student.aau.dk/sheide17/:/home/nanodisco/dataset  nd_env # GPU server

######################################################################################################
# NanoDisco current difference
######################################################################################################
ND_difference() {
  difference=$(dirname $1)/difference
  echo "PRE    $1"
  echo "DIF    $difference"
  echo "-w     $3"
  echo "-n     $4"
  echo "REF    $2"
  
  if [ ! -s $difference ]; then
    nanodisco difference -nj 50 -nc 2 -p 10 \
    -i $1 \
    -o $difference \
    -w $3 \
    -n $4 \
    -r $2
  fi
}
# BIO-server path "/shared-nfs/SH/samples/"
# AI-cloud path "projects/modification_binning/"

ND_difference "$nd_path/shared-nfs/SH/samples/zymo/nanodisco/preprocess" \
  "$nd_path/shared-nfs/SH/samples/zymo/nanodisco/assembly_100x.fasta" \
  zymo_wga \
  zymo_nat

nanodisco merge -d $nd_path/shared-nfs/SH/samples/zymo/nanodisco/difference \
  -o $nd_path/shared-nfs/SH/samples/zymo/nanodisco/difference_merge \
  -b zymo


nanodisco coverage -b $nd_path/shared-nfs/SH/samples/zymo/nanodisco/preprocess/zymo_wga.sorted.bam \
  -r $nd_path/shared-nfs/SH/samples/zymo/nanodisco/assembly_100x.fasta \
  -o $nd_path/shared-nfs/SH/samples/zymo/nanodisco/coverage
nanodisco coverage -b $nd_path/shared-nfs/SH/samples/zymo/nanodisco/preprocess/zymo_nat.sorted.bam \
  -r $nd_path/shared-nfs/SH/samples/zymo/nanodisco/assembly_100x.fasta \
  -o $nd_path/shared-nfs/SH/samples/zymo/nanodisco/coverage


nanodisco profile -p 100 \
  -r $nd_path/shared-nfs/SH/samples/zymo/nanodisco/assembly_100x.fasta \
  -d $nd_path/shared-nfs/SH/samples/zymo/nanodisco/difference_merge/zymo_difference.RDS \
  -w $nd_path/shared-nfs/SH/samples/zymo/nanodisco/coverage/zymo_wga.cov \
  -n $nd_path/shared-nfs/SH/samples/zymo/nanodisco/coverage/zymo_nat.cov \
  -b zymo \
  -a "all" \
  -o $nd_path/shared-nfs/SH/samples/zymo/nanodisco/profile \
  --min_contig_len 10000 \
  -c 5
# if h5dump -pH "/shared-nfs/SH/samples/zymo/NAT/fast5_multi_100x_filt/fast5_0.fast5" | grep -q vbz; then
  # echo "It is VBZ-compressed"
# else
  # echo "It is not VBZ-compressed"
# fi





######################################################################################################
# NanoDisco AD sample
######################################################################################################

nanodisco_pipeline.sh -n "/shared-nfs/SH/samples/AD/NAT/fast5_pooled_reads_filt" \
  -w "/shared-nfs/SH/samples/AD/WGA/fast5_pooled_reads_filt" -a "/shared-nfs/SH/samples/AD/polished_assembly.fasta" \
  -o "/shared-nfs/SH/samples/AD/nanodisco" \
  -t 2
