#!/bin/bash

######################################################################################################
# Initialisation
######################################################################################################

# Modules


# Help function for inputs
helpFunction () {
  echo ""
  echo "Basic usage: bash $0 -p PROFILE -a ASSEMBLY -o OUT"
  echo "  -p|--profile     methylation profile generated by Nanodisco" 
  echo "  -a|--assembly    assembly with contigs for binning"
  echo "  -o|--out         output directionary"
  echo "  --embedding      method for feature embedding (Autoencoder)"
  echo "  --threads        number of threads to use (default:10)"
  
  exit 1 # Exit script after printing help
}

# Optional deafaults parameter settings
threads=10

# Get parameters
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      helpFunction
      shift ;;
    -p|--profile)
      profile=$2
      shift 2 ;;
    -a|--assembly)
      assembly=$2
      shift 2 ;;
    -o|--out)
      out=$2
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
    --embedding)
      embedding=$2
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
if [ -z "$profile" ] || [ -z "$assembly" ]; then
  helpFunction
fi

# Defineing log file and function
mkdir -p $out
log () {
  echo $(date +%d-%m-%Y%t%k:%M:%S) " $1"
  echo $(date +%d-%m-%Y%t%k:%M:%S) " $1" >> $out/log_meth_bin.txt
}
echo " " >> $out/log_meth_bin.txt




######################################################################################################
# Prepare methylation features
######################################################################################################

mes="  Feature selection"
if [[ ! -s "$out/features.tsv" ]]; then
  log "$mes"
  
  feature_preprocessing.R -p $profile -a $assembly -o $out

  else log "SKIP: $mes"
fi
features=$out/features.tsv

######################################################################################################
# Cluster contigs based on methylation features
######################################################################################################

mes="  Embedding"
if [[ ! -s "$out/AE/binnning_AE_representation.tsv" ]]; then
  log "$mes"
  
  conda activate py38_pyto
  python3 AE_modbin.py $features $out/AE

  else log "SKIP: $mes"
fi



mes="  Clustering"
if [[ ! -s "$out/bins_AE" ]]; then
  log "$mes"
  
  clustering.R -e $out/AE/binnning_AE_representation.tsv -a $assembly -o $out
  
  else log "SKIP: $mes"
fi



mes="  Checking bin quality"
if [[ ! -s "$out/bins_AE/summary.txt" ]]; then
  log "$mes"
  
  bin_QC.sh -i "$out/bins_AE" \
    -o "$out/bins_AE" \
    -q $fastq \
    -x ".fasta"
  
  else log "SKIP: $mes"
fi
