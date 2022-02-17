#!/bin/bash

######################################################################################################
# Initialisation
######################################################################################################
# Stop if execution of command fails

# Modules
mod_parallel=parallel/20190122-foss-2018a                 # https://www.gnu.org/software/bash/manual/html_node/GNU-Parallel.html
mod_nanoplot=NanoPlot/1.38.0-foss-2020b                   # https://github.com/wdecoster/NanoPlot
mod_nanofilt=nanofilt/2.6.0-foss-2020b-Python-3.8.6       # https://github.com/wdecoster/nanofilt

# Help function for inputs
helpFunction () {
  echo ""
  echo "Basic usage: bash $0 -f Fast5 -q Fastq"
  echo "  -q | --fastq     name of fastq file to write into (if it exist the filtration and QC is only performed)"
  echo "  -t | --threads   number of threads to use (default:10)"
  echo "  --filt           include flag to perform filtration (othervise, only quast QC)"
  echo "  --qc_only         Only performs QC and filtration on specified fastq file"
  echo "  --min_qual       minimum mean quality of read if filtering (default:5)"
  echo "  --min_length     minimum mean length of read if filtering (default:1000)"
  
  exit 1 # Exit script after printing help
}

# Optional deafaults parameter settings
threads=10
min_qual=5
min_length=1000
filt=0
qc_only=0
# Get parameters
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      helpFunction
      shift ;;
    -f|--fast5)
      fast5=$2
      shift 2 ;;
    -q|--fastq)
      fastq=$2
      shift 2 ;;
    -t|--threads)
      threads=$2
      shift 2 ;;
    --filt)
      filt=1
      shift ;;
    --min_length)
      min_length=$2
      shift 2 ;;
    --min_qual)
      min_qual=$2
      shift 2 ;;
    --qc_only)
      qc_only = 1
      shift ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1 ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift ;;
  esac
done
eval set -- "$PARAMS"

# Required parameter check
shift "$(( OPTIND - 1 ))"
if [ -z "$fastq" ]; then
  helpFunction
fi

# Defineing log file and function
log () {
  echo $(date +%d-%m-%Y%t%k:%M:%S) " $1" >> $(dirname $fastq)/log_prep.txt
}
echo " " >> $(dirname $fastq)/log_prep.txt

######################################################################################################
# Extracting fastq
######################################################################################################
# The fastq are also availible in the AAU_CS folderm but i didn't have reading right, so i just \n
# extracted from the fast5 reads 
log "Starting script"


#single_fastq_dir=${fastq%.*}
#mkdir -p $single_fastq_dir
#module purge
#extract_fastq () {
#    log "[INIT] Extracting fastq from $fast5 to $single_fastq_dir"
#    python fast5_to_fastq.py $fast5 "$single_fastq_dir" --jobs 400 --batch_size 4000
#    log "[DONE]"
#}


#extract_fastq \
#	$fast5 \
#	$fastq

######################################################################################################
# Quality check of fastq files
######################################################################################################

nanoplot () {
  dir=$(dirname $1)
  nanoplot_dir=$dir/$(basename $1 .fastq)_nanoplot
  if [[ ! -s "$nanoplot_dir/NanoPlot-report.html" ]]; then
    log "[INIT] Quality check on $1"
    mkdir -p "$nanoplot_dir"
    module load $mod_nanoplot
    NanoPlot --raw --tsv_stats -t $threads --fastq $fastq -o $nanoplot_dir
    module purge
    log "[DONE]"
    else log "[SKIP] Quality check already performed on $1 (see $nanoplot_dir)"
  fi
}

nanoplot \
  $fastq

######################################################################################################
# Filtering and quality check of fastq
######################################################################################################

fastq_filt=${fastq%.*}_filt.fastq
if [[ $filt == 1 ]] && [[ ! -s $fastq_filt ]]; then
  log "[INIT] Filtering fastq with min_qual: $min_qual and min_length: $min_length"
  module load $mod_nanofilt
  NanoFilt -q $min_qual -l $min_length $fastq > $fastq_filt
  log "[DONE]"
  else log "[SKIP] Filtering already performed"
fi

nanoplot \
  $fastq_filt







