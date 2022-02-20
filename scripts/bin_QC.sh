#!/bin/bash

######################################################################################################
# Initialisation
######################################################################################################
# Stop if execution of command fails
# Modules
mod_checkm=CheckM/1.1.2-foss-2018a-Python-3.6.4         # https://github.com/Ecogenomics/CheckM/wiki
mod_samtools=SAMtools/1.11-foss-2020b                   # http://www.htslib.org/doc/
mod_minimap=Minimap2/2.17-foss-2020b                    
mod_gtdbtk=GTDBTk/1.5.0-foss-2020b-Python-3.8.6

# Optional deafaults parameter settings
threads=10
extension="fasta"
coverage=0

# Help function for inputs
helpFunction () {
  echo "
Basic usage: bash $0 -i [Input] -o [Output] 
  -i | --bins        path to bins with fasta sequences
  -o | --out         output directionary
  -q | --fastq       path to fastq folders, which should be named after basename of fasta file, e.g. as [q]/basename/basename.fastq
  -x | --extension   extension of fasta files (deafault: $extension)
  -t | --threads     number of threads to use (default: $threads)
  -c | --coverage    include flag to caluclate coverage (require --fastq to be specified)
"
  exit 1 # Exit script after printing help
}

# Get parameters
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      helpFunction
      shift ;;
    -i|--bins)
      bins=$2
      shift 2 ;;
    -o|--out)
      out=$2
      shift 2 ;;
    -x|--extension)
      extension=$2
      shift 2 ;;
    -t|--threads)
      threads=$2
      shift 2 ;;
    -q|-fastq)
      fastq=$2
      shift 2 ;;
    -c|--coverage)
      coverage=1
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
if [[ -z "$bins" ]] || [[ -z "$out" ]]; then
  helpFunction
fi

# Defineing log file and function
log () {
  echo $(date +%d-%m-%Y%t%k:%M:%S) " $1" >> $(dirname $bins)/log_bin_check.txt
}

######################################################################################################
# CheckM
######################################################################################################
module purge
mkdir $out

# CheckM
if [[ ! -s $out/summary.txt ]]; then
  log "[STARTED] CheckM"
  module load $mod_checkm
  mkdir $bins/.temp
  checkm lineage_wf $bins $out -x $extension -t $threads --tmpdir $bins/.temp --file summary.txt
  checkm qa $out/lineage.ms $out > $out/summary.txt
  rm -r $bins/.temp
  module purge
  log "[DONE] CheckM"
  else log "[SKIP] CheckM already performed on bins (see $out)"
fi

######################################################################################################
# GTDBTk
######################################################################################################

if [[ ! -s $out/GTDBTk ]]; then
  log "[STARTED] GTDBTk"
  module load $mod_gtdbtk
  gtdbtk classify_wf  -x $extension \
    --cpus $threads \
    --genome_dir $bins \
    --out_dir $out/GTDBTk
  module purge
  log "[DONE] CheckM"
  else log "[SKIP] CheckM already performed on bins (see $out)"
fi


######################################################################################################
# Coverage
######################################################################################################

