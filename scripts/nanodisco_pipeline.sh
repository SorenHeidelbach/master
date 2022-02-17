#!/bin/bash

######################################################################################################
# Initialisation
######################################################################################################

# Modules

# Optional deafaults parameter settings
threads=10

# Help function for inputs
helpFunction () {
  echo ""
  echo "Basic usage: bash $0 -n NAT.fast5 -w WGA.fast5 -a ASSEMBLY.fasta"
  echo "  -n | --fast5-nat  NAT fast5 folder"
  echo "  -w | --fast5-wga  WGA fast5 folder"
  echo "  -a | --assembly   assembly with contigs"
  echo "  -o | --out        output folder for nanodisco analysis"
  echo "  -t | --threads    number of threads to use"
  
  
  exit 1 # Exit script after printing help
}

# Get parameters
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      helpFunction
      shift ;;
    -n|--fast5-nat)
      export fast5_nat=$2
      shift 2 ;;
    -w|--fast5-wga)
      export fast5_wga=$2
      shift 2 ;;
    -a|--assembly)
      export assembly=$2
      shift 2 ;;
    -o|--out)
      export wd=$2
      shift 2 ;;
    -t|--threads)
      export threads=$2
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

# Required parameter check
shift "$(( OPTIND - 1 ))"
if [ -z "$fast5_nat" ] || [ -z "$fast5_wga" ] || [ -z "$assembly" ]; then
  helpFunction
fi
if [ -z "$wd" ]; then
  export wd=$(dirname $fast5)
fi

# Defineing log file and function
mkdir -p $wd
log () {
  if [[ $2 == 0 ]]; then
    echo ""
  else
    echo $1
    echo $(date +%d-%m-%Y%t%k:%M:%S) " $1" >> $wd/log.txt
  fi
}
export -f log
echo " " >> $wd/log.txt
# =================================================================================================
# =================================================================================================

# Nanodisco require fast5 to be gzip compressed and not VBZ-compressed
conda activate py38
compress_fast5_fun() {
  compress_fast5 -t $threads --recursive \
    -i $1 \
    -s "$1"_gzip \
    -c gzip
}

mes="Compress NAT fast5 to gzip"
if [[ ! -s ${fast5_nat}_gzip ]]; then
  log "$mes"
  
  compress_fast5_fun $fast5_nat
  
  else log "SKIP: $mes"
fi
export fast5_nat_gzip=${fast5_nat}_gzip

mes="Compress WGA fast5 to gzip"
if [[ ! -s ${fast5_wga}_gzip ]]; then
  log "$mes"
  
  compress_fast5_fun $fast5_wga
  
  else log "SKIP: $mes"
fi
export fast5_wga_gzip=${fast5_wga}_gzip

# Initiate Nanodisco container
singularity exec --no-home -B /:/home/nanodisco/dataset /shared-nfs/SH/software/nanodisco/nd_env bash /home/nanodisco/dataset/shared-nfs/SH/code/scripts/nanodisco_commands.sh
# wd=/home/nanodisco/dataset/$wd

#Preprocessing
# mes="Preprocess WGA "
# if [[ ! -s ${fast5_wga}_gzip ]]; then
  # log "$mes"
  
#  WGA
  # nanodisco preprocess -p $threads \
    # -f $fast5_nat_gzip \
    # -o $wd/preprocess \
    # -r $assembly \
    # -s nat
  
  # else log "SKIP: $mes"
# fi

# mes="Preprocess NAT "
# if [[ ! -s ${fast5_wga}_gzip ]]; then
  # log "$mes"
  
#  NAT
  # nanodisco preprocess -p $threads \
    # -f $fast5_wga_gzip   \
    # -o $wd/preprocess \
    # -r $assembly \
    # -s wga
  
  # else log "SKIP: $mes"
# fi

# Current difference
# mes="Current differece"
# if [[ ! -s ${fast5_wga}_gzip ]]; then
  # log "$mes"
  
  # nanodisco difference -nj 50 -nc 2 -p 10 \
    # -i $wd/preprocess \
    # -o $wd/difference \
    # -w wga \
    # -n nat \
    # -r $assembly
  
  # else log "SKIP: $mes"
# fi

# nanodisco merge -d $wd/difference \
  # -o $wd/difference_merge \
  # -b sample

# Coverage
# mes="Coverage WGA"
# if [[ ! -s $wd/coverage/wga.cov ]]; then
  # log "$mes"
  
  # nanodisco coverage -b $wd/preprocess/zymo_wga.sorted.bam \
    # -r $assembly \
    # -o $wd/coverage
  
  # else log "SKIP: $mes"
# fi
# mes="Coverage NAT"
# if [[ ! -s $wd/coverage/nat.cov ]]; then
  # log "$mes"
  
  # nanodisco coverage -b $wd/preprocess/zymo_nat.sorted.bam \
    # -r $assembly \
    # -o $wd/coverage
  
  # else log "SKIP: $mes"
# fi

# Methylation profile
# mes="Coverage NAT"
# if [[ ! -s $wd/profile/methylation_profile_sample.RDS]]; then
  # log "$mes"
  
  # nanodisco profile -p $threads \
    # -r $assembly \
    # -d $wd/difference_merge/zymo_difference.RDS \
    # -w $wd/coverage/zymo_wga.cov \
    # -n $wd/coverage/zymo_nat.cov \
    # -b sample \
    # -a "all" \
    # -o $wd/profile \
    # --min_contig_len 5000 \
    # -c 5
  
  # else log "SKIP: $mes"
# fi

