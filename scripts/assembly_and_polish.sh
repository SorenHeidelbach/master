#!/bin/bash

######################################################################################################
# Initialisation
######################################################################################################
# Stop if execution of command fails

# Modules
mod_racon=Racon/1.3.3-pikachu-foss-2018a                  # https://denbi-nanopore-training-course.readthedocs.io/en/latest/polishing/medaka/racon.html
mod_medaka=Medaka/1.2.3-foss-2020b                        # https://github.com/nanoporetech/medaka
mod_minimap=Minimap2/2.17-foss-2020b                      # https://github.com/lh3/minimap2
mod_flye=Flye/2.9-GCC-10.2.0                              # https://github.com/fenderglass/Flye
mod_quast=QUAST/4.6.3-foss-2018a-Python-3.6.4             # http://quast.sourceforge.net/docs/manual.html#sec2

# Optional deafaults parameter settings
threads=10
genome_size=50m

# Help function for inputs
helpFunction () {
  echo ""
  echo "Basic usage: bash $0 -q Fastq"
  echo "  -q | --fastq     name of fastq file with all reads (e.g. /home/analysis/sample/sample1.fastq)"
  echo "  -o | --out-dir   directionary to output assembly and polished assemblies"
  echo "  -t | --threads   number of threads to use (default:$threads)"
  echo "  --genome-size    passed directly to flye (default:$genome_size)"
  
  exit 1 # Exit script after printing help
}

# Get parameters
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      helpFunction
      shift ;;
    -q|--fastq)
      fastq=$2
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

# Required parameter check
shift "$(( OPTIND - 1 ))"
if [ -z "$fastq" ]; then
  helpFunction
fi
if [ -z "$out_dir" ]; then
  out_dir=$(dirnmae $fastq)/assembly
fi

# Defineing log file and function
mkdir -p $out_dir
log () {
  echo $(date +%d-%m-%Y%t%k:%M:%S) " $1" >> $out_dir/log_assembly.txt
}
echo " " >> $out_dir/log_assembly.txt

######################################################################################################
# Assembly
######################################################################################################
log "Start of script"
module purge

# Flye assembly
if [ ! -s $out_dir/flye/assembly.fasta ]; then
  module load $mod_flye
  log "[INIT] Starting assembly, genome-size:$genome_size"
  mkdir -p $out_dir/flye
  flye --nano-raw $fastq --genome-size $genome_size --out-dir $out_dir/flye --meta --threads $threads
  log "[DONE]"
  module purge
  else log "[SKIP] Assembly already present in $out_dir/flye"
fi
assembly=$out_dir/flye/assembly.fasta


# Quast
if [[ ! -s $out_dir/quast ]]; then
  log "[INIT] Flye QC with quast"
  mkdir $out_dir/quast
  module load $mod_quast
  if [[ ! -s $out_dir/quast/alignment.sam ]]; then
    log "[INIT] Mapping $fastq to $assembly"
    module load $mod_minimap
    minimap2 -ax map-ont -t $threads \
      $assembly \
      $fastq > \
      $out_dir/quast/alignment.sam
    log "[DONE]"
    else log "[SKIP] Alignment to reference already present: $out_dir/quast/alignment.sam"
  fi
  quast.py $assembly \
    --sam $out_dir/quast/alignment.sam \
    -o $out_dir/quast
  module purge
  log "[DONE]"
  else log "[SKIP] Flye QC, Quast folder already detected"
fi



######################################################################################################
# Polishing
######################################################################################################

# Racon
run_racon () {
  if [ ! -s $out_dir/racon$2/assembly.fasta ]; then
  log "[INIT] Racon polish on $1"
  module load $mod_racon
  module load $mod_minimap
  mkdir -p $out_dir/racon$2
  minimap2 -x map-ont -t $threads \
    $1 \
    $fastq > \
    $out_dir/racon$2/alignment.paf
  racon -t $threads \
    $fastq \
    $out_dir/racon$2/alignment.paf \
    $1 > \
    $out_dir/racon$2/assembly.fasta 
  module purge
  log "[DONE] Polished assembly in $out_dir/racon$2"
  else log "[SKIP] $o_dir/racon$2/assembly.fasta was identified and therefore Racon polish nr. $2 skipped"
  fi
}

run_racon $assembly 1
run_racon $out_dir/racon1/assembly.fasta 2
run_racon $out_dir/racon2/assembly.fasta 3
run_racon $out_dir/racon3/assembly.fasta 4

# Medaka
run_medaka () {
  if [ ! -s $out_dir/medaka$3/consensus.fasta ]; then
  log "[INIT] Polishing with medaka on $2"
  module load $mod_medaka
  medaka_consensus -i $2 -d $1 -o $out_dir/medaka$3 -t $threads
  log "[DONE]"
  else log "[SKIP] Medaka ($out_dir/medaka$3/consensus.fasta present)"
  fi
}
run_medaka $fastq $out_dir/racon4/assembly.fasta "_1"
run_medaka $fastq $out_dir/medaka1/consensus.fasta "_2"


######################################################################################################
# Assembly QC
######################################################################################################










