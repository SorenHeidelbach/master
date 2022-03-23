#!/bin/bash

######################################################################################################
# Initialisation
######################################################################################################
# Stop if execution of command fails


# Modules
mod_cutadapt=cutadapt/3.4-GCCcore-10.2.0-Python-3.8.6
mod_minimap=Minimap2/2.17-foss-2020b                      # https://github.com/lh3/minimap2
mod_nanofilt=nanofilt/2.6.0-foss-2020b-Python-3.8.6
mod_nanoplot=NanoPlot/1.38.0-foss-2020b                   # https://github.com/wdecoster/NanoPlot


# Optional deafaults parameter settings
threads=10
genome_size=50m
FAST5_input=1
# Help function for inputs
helpFunction () {
  echo ""
  echo "Basic usage: bash $0 -o OUTPUT_PATH --fastq FASTQ"
  echo "  --fastq             Name of directionary with basecalled fastq"
  echo "  -t | --threads      Number of threads to use (default:$threads)"
  echo "  -o | --out          Output directionary"
  
  exit 1 # Exit script after printing help
}

# Get parameters
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      helpFunction
      shift ;;
    --fastq)
      fastq=$2
      shift 2 ;;
    -o|--out-dir)
      out_dir=$2
      shift 2 ;;
    -t|--threads)
      threads=$2
      shift 2 ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1 ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift ;;
  esac
done
eval set -- "$PARAMS"

# Defineing log file and function
mkdir -p $out_dir
log () {
  echo $(date +%d-%m-%Y%t%k:%M:%S) " $1" >> $out_dir/log.txt
}
echo " " >> $out_dir/log.txt



######################################################################################################
# Preprocessing of reads
######################################################################################################
#-- QC of extracted raw reads --
fun_nanoplot() {
echo IN $1 - OUT $(dirname $1)/nanoplot_$(basename ${1%.*})
  module load $mod_nanoplot
  NanoPlot --raw --tsv_stats -t $threads --fastq $1 -o $(dirname $1)/nanoplot_$(basename ${1%.*})
  module purge
}

OUT=$(dirname $fastq)/nanoplot_$(basename ${fastq%.*})
if [ ! -s $OUT ]; then
  log "[init] Nanoplot on $fastq"
  fun_nanoplot $fastq
  log "[done] "
  else log "[skip] $OUT already exists"
fi


#-- Removal of adapters --
fun_cutadapt () {
  module load $mod_cutadapt
  echo IN"  "$1
  echo OUT" "${1%.*}_trim.fastq
  cutadapt -g "XGTGAATGTACTTCGTTCAGTTACGTATTGCT;min_overlap=10" -a "GCAATACGTAACTGAACGAAGTX" -e 0.2 -j 20 -o ${1%.*}_trim.fastq $1 > $(dirname $1)/trimming_summary.txt
  module purge
}

OUT=${fastq%.*}_trim.fastq
if [ ! -s $OUT ]; then
  log "[init] Removing adapters from $fastq"
  fun_cutadapt $fastq
  log "[done] Saved to $OUT"
  else log "[skip] $OUT exists"
fi
fastq_trim=$OUT


#-- Filtration of reads --
fun_nanofilt () {
  fastq_filt=${1%.*}_filt.fastq
  module load $mod_nanofilt
  NanoFilt -q 10 -l 1000 $1 > $fastq_filt
  module purge
  fun_nanoplot $fastq_filt
}

OUT=${fastq_trim%.*}_filt.fastq
if [ ! -s $OUT ]; then
  log "[init] Filtration of $fastq_trim"
  fun_nanofilt $fastq_trim
  log "[done] "
  else log "[skip] $OUT exists"
fi
fastq_filt=$OUT

#-- QC of filtered and trimmed reads --
OUT=$(dirname $fastq)/nanoplot_$(basename ${fastq_filt%.*})
if [ ! -s $OUT ]; then
  log "[init] Nanoplot on $fastq_filt"
  fun_nanoplot $fastq_filt
  log "[done] "
  else log "[skip] $OUT already exists"
fi

