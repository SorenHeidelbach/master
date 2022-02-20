#!/bin/bash

######################################################################################################
# Initialisation
######################################################################################################
# Stop if execution of command fails
set -e pipefail

# Modules
mod_racon=Racon/1.3.3-pikachu-foss-2018a                  # https://denbi-nanopore-training-course.readthedocs.io/en/latest/polishing/medaka/racon.html
mod_medaka=Medaka/1.2.3-foss-2020b                        # https://github.com/nanoporetech/medaka
mod_minimap=Minimap2/2.17-foss-2020b                      # https://github.com/lh3/minimap2
mod_nanoplot=NanoPlot/1.38.0-foss-2020b                   # https://github.com/wdecoster/NanoPlot
mod_flye=Flye/2.9-GCC-10.2.0                              # https://github.com/fenderglass/Flye
mod_quast=QUAST/4.6.3-foss-2018a-Python-3.6.4             # http://quast.sourceforge.net/docs/manual.html#sec2

# Optional deafaults parameter settings
threads=10
genome_size=50m
FAST5_input=1
# Help function for inputs
helpFunction () {
  echo ""
  echo "Basic usage: bash $0 -o OUTPUT_PATH -w WGA_FAST5 -n NAT_FAST5 OR -o OUTPUT_PATH --WGA_fastq WGA_FASTQ --NAT_fastq NAT_FASTQ"
  echo "  -w | --WGA_fast5    Name of directionary with WGA basecalled fast5"
  echo "  -n | --NAT_fast5    Name of directionary with NAT basecalled fast5"
  echo "  --WGA_fastq         Name of directionary with NAT basecalled fastq"
  echo "  --NAT_fastq         Name of directionary with NAT basecalled fastq"
  echo
  echo "  -t | --threads      Number of threads to use (default:$threads)"
  echo "  -o | --out          Output directionary"
  echo "  --genome-size       Passed directly to flye (default:$genome_size)"
  
  exit 1 # Exit script after printing help
}

# Get parameters
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      helpFunction
      shift ;;
    -w|--WGA_fast5)
      WGA_fast5=$2
      shift 2 ;;
    -n|--NAT_fast5)
      NAT_fast5=$2
      shift 2 ;;
    --WGA_fastq)
      WGA_fastq=$2
      shift 2 ;;
    --NAT_fastq)
      NAT_fastq=$2
      shift 2 ;;
    -o|--out-dir)
      out_dir=$2
      shift 2 ;;
    -t|--threads)
      threads=$2
      shift 2 ;;
    --genome-size)
      genome_size=$2
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


# Check if fast5 or fastq are specified
shift "$(( OPTIND - 1 ))"
if [ -z "$WGA_fast5" ] && [ -z "$NAT_fast5" ]; then
  FAST5_input=0
  if [ -z "$WGA_fastq" ] && [ -z "$NAT_fastq" ]; then
    helpFunction
  fi
fi

# Check for output folder
if [ -z $out_dir ]; then
  helpFunction
fi

# Defineing log file and function
mkdir -p $out_dir
log () {
  echo $(date +%d-%m-%Y%t%k:%M:%S) " $1" >> $out_dir/log_assembly.txt
}
echo " " >> $out_dir/log.txt


######################################################################################################
# Extract Fastq from fast5
######################################################################################################


if [ $FAST5_input == 1 ]; then
  log "Fast5 are checked for basecalled reads"
  # WGA extraction
  mkdir -p "$out_dir/WGA"
  OUT="$out_dir/WGA/reads.fastq"
  if [ ! -s $OUT ]; then
    log "[init] Fastq are extracted from $WGA_FAST5"
    python /shared-nfs/SH/code/scripts/fast5_to_fastq.py $WGA_FAST5 $OUT --jobs $threads --batch_size 4000
    log "[done] Fastqs are in $OUT"
    else log "[skip] Fastq file is already present ($OUT)"
  fi
  WGA_fastq="$out_dir/WGA/reads.fastq"
  # NAT fastq extraction
  mkdir -p "$out_dir/NAT"
  OUT="$out_dir/NAT/reads.fastq"
  if [ ! -s $OUT ]; then
    log "[init] Fastq are extracted from $NAT_FAST5"
    python /shared-nfs/SH/code/scripts/fast5_to_fastq.py $WGA_FAST5 $OUT --jobs $threads --batch_size 4000
    log "[done] Fastqs are in $OUT"
    else log "[skip] Fastq file is already present ($OUT)"
  fi
  NAT_fastq="$out_dir/NAT/reads.fastq"
fi


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

# WGA
OUT=$(dirname $WGA_fastq)/nanoplot_$(basename ${WGA_fastq%.*})
if [ ! -s $OUT ]; then
  log "[init] Nanoplot on $WGA_fastq"
  fun_nanoplot $WGA_fastq
  log "[done] "
  else log "[skip] $OUT already exists"
fi
# NAT
OUT=$(dirname $NAT_fastq)/nanoplot_$(basename ${NAT_fastq%.*})
if [ ! -s $OUT ]; then
  log "[init] Nanoplot on $NAT_fastq"
  fun_nanoplot $NAT_fastq
  log "[done] "
  else log "[skip] $OUT already exists"
fi



#-- Removal of adapters --
fun_cutadapt () {
  module load cutadapt/3.4-GCCcore-10.2.0-Python-3.8.6
  echo IN"  "$1
  echo OUT" "${1%.*}_trim.fastq
  cutadapt -g "XGTGAATGTACTTCGTTCAGTTACGTATTGCT;min_overlap=10" -a "GCAATACGTAACTGAACGAAGTX" -e 0.2 -j 20 -o ${1%.*}_trim.fastq $1 > $(dirname $1)/trimming_summary.txt
  module purge
}

#WGA
OUT=${WGA_fastq%.*}_trim.fastq
if [ ! -s $OUT ]; then
  log "[init] Removing adapters from $WGA_fastq"
  fun_cutadapt $WGA_fastq
  log "[done] Saved to $OUT"
  else log "[skip] $OUT exists"
fi
WGA_fastq_trim=$OUT

#NAT
OUT=${NAT_fastq%.*}_trim.fastq
if [ ! -s $OUT ]; then
  log "[init] Removing adapters from $NAT_fastq"
  fun_cutadapt $NAT_fastq
  log "[done] Saved to $OUT"
  else log "[skip] $OUT exists"
fi
NAT_fastq_trim=$OUT



#-- Filtration of reads --
fun_nanofilt () {
  fastq_filt=${1%.*}_filt.fastq
  module load nanofilt/2.6.0-foss-2020b-Python-3.8.6
  NanoFilt -q 7 -l 1000 $1 > $fastq_filt
  module purge
  fun_nanoplot $fastq_filt
}

# WGA 
OUT=${WGA_fastq_trim%.*}_filt.fastq
if [ ! -s $OUT ]; then
  log "[init] Filtration of $WGA_fastq_trim"
  fun_nanofilt $WGA_fastq_trim
  log "[done] "
  else log "[skip] $OUT exists"
fi
WGA_fastq_filt=$OUT

# NAT
OUT=${NAT_fastq_trim%.*}_filt.fastq
if [ ! -s $OUT ]; then
  log "[init] Filtration of $NAT_fastq_trim"
  fun_nanofilt $NAT_fastq_trim
  log "[done] "
  else log "[skip] $OUT exists"
fi
NAT_fastq_filt=$OUT







