#!/bin/bash

######################################################################################################
# Initialisation
######################################################################################################

# Modules
mod_minimap=Minimap2/2.17-foss-2020b                      # https://github.com/lh3/minimap2
mod_metabat=MetaBAT/2.12.1-foss-2018a
mod_maxbin=MaxBin/2.2.7-foss-2018a-Perl-5.26.1
mod_samtools=SAMtools/1.11-foss-2020b

# Optional deafaults parameter settings
threads=10

# Help function for inputs
helpFunction () {
  echo ""
  echo "Basic usage: bash $0 -a ASSEMBLY.fasta -q READS.fastq"
  echo "  -a | --assembly  assembly with contigs"
  echo "  -q | --fastq     name of fastq file with all reads"
  echo "  -o | --out-dir   directionary to output binning results (deafault:same level as assembly)"
  echo "  -t | --threads   number of threads to use (default:$threads)"
  
  exit 1 # Exit script after printing help
}

# Get parameters
PARAMS=""
while (( "$#" )); do
  case "$1" in
    -h|--help)
      helpFunction
      shift ;;
    -a|--assembly)
      metagenome=$2
      shift 2 ;;
    -q|--fastq)
      fastq=$2
      shift 2 ;;
    -o|--out-dir)
      wd=$2
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

# Required parameter check
shift "$(( OPTIND - 1 ))"
if [ -z "$fastq" ] || [ -z "$metagenome" ]; then
  helpFunction
fi
if [ -z "$wd" ]; then
  wd=$(dirname $metagenome)
fi

# Defineing log file and function
mkdir -p $wd
log () {
  if [[ $2 == 0 ]]; then
    echo ""
  else
    echo $1
    echo $(date +%d-%m-%Y%t%k:%M:%S) " $1" >> $wd/log_binning.txt
  fi
}
echo " " >> $wd/log_binning.txt

bam=$wd/$(basename ${fastq%.fa*})_to_$(basename ${metagenome%.fa*}).bam

# =================================================================================================
# =================================================================================================

# =================================================================================================
# Read Mapping

mes="Read mapping"
if [[ ! -s $bam ]]; then
  log "$mes"
  module purge
  module load $mod_minimap
  module load $mod_samtools
  minimap2 -ax map-ont -t $threads $metagenome \
    $fastq | \
    samtools view --threads $threads -Sb | samtools sort --threads $threads -o $bam
  module purge
  
  else log "SKIP: $mes"
fi
# =================================================================================================
# Coverage files

calculate_coverage_from_bam() {
  log "Coverage extraction on $1"
  module purge
  module load $mod_samtools
  
  mes="  Calculating .depth file from $bam"
  if [[ ! -s ${1%.bam}.depth ]]; then
    log "$mes"
    # command
    samtools depth -a $1 > ${1%.bam}.depth
    #
    else log "SKIP: $mes"
  fi
  
  mes="  Extraction mean contig coverage"
  if [[ ! -s ${1%.bam}.cov.csv ]]; then
    log "$mes"
    # command
    echo "scaffold,coverage" > ${1%.bam}.cov.csv
    awk -F '\t' '{a[$1] += $3;I[$1]++} END{for (i in a) print i,",",a[i]/I[i]}' ${1%.bam}.depth >> ${1%.bam}.cov.csv
    sed -i "s/ , /,/"  ${1%.bam}.cov.csv
    #
    else log "SKIP: $mes"
  fi  
  
  mes="  Converting to tsv"
  if [[ ! -s ${1%.bam}.cov.tsv ]]; then
    log "$mes"
    # command
    tail -n+2 ${bam%.bam}.cov.csv | sed "s/,/\t/" > ${1%.bam}.cov.tsv
    #
    else log "SKIP: $mes"
  fi
  module purge
}
calculate_coverage_from_bam $bam

# Subset to contigs that have at least one read mapping (required for maxbin)
awk '{
if (NR == FNR) {a[">"$1];next}
if($1 ~ />/) {
  if ($1 in a) write_out=1
  else write_out=0
  }
if (write_out == 1) print $1}' ${bam%.bam}.cov.tsv $metagenome > ${metagenome%.fasta}_mapped.fasta

# =================================================================================================
# Checkm summary file conversion to easily readable file
checkm_summary_to_tsv() {
  awk '!/[\[\-]/ {print$0}' $1 | \
    sed -e "s:\#\s::g" | \
    sed "s/Bin ID/bin_id/" | \
    sed "s/Marker lineage/marker_lineage/" | \ 
    sed "s/Strain heterogeneity/strain_heterogeneity/" > $(dirname $1)/summary.tsv
} 

# =================================================================================================
# Metabat2
log "Binning"
mes="  Metabat2"
if [[ ! -s $wd/metabat2/metabat.1.fa ]]; then
  log "$mes"
  
  module purge
  module load $mod_metabat
  metabat2 -i $metagenome $bam -o $wd/metabat2/metabat
  module purge
  
  else log "SKIP: $mes"
fi
# Check bin quality
mes="  Metabat checkm"
if [[ ! -s $wd/metabat2/summary.tsv ]]; then
  log "$mes"
  
  bash /shared-nfs/SH/code/scripts/cluster_quality_check.sh -i "$wd/metabat2" \
    -o "$wd/metabat2" \
    -q $fastq \
    -x ".fa"
  checkm_summary_to_tsv "$wd/metabat2/summary.txt"
  
  else "SKIP: $mes"
fi
# =================================================================================================
# Maxbin
mes="  Maxbin"
if [[ ! -s $wd/maxbin/maxbin.001.fasta ]]; then
  log "$mes"
  
  module purge
  module load $mod_maxbin
  mkdir -p $wd/maxbin
  run_MaxBin.pl -contig ${metagenome%.fasta}_mapped.fasta -abund ${bam%.bam}.cov.tsv -out $wd/maxbin/maxbin -thread $threads
  module purge
  
  else "SKIP: $mes"
fi
# Check bin quality
mes="  Maxbin checkm"
if [[ ! -s $wd/maxbin/summary.tsv ]]; then
  log "$mes"
  
  bash /shared-nfs/SH/code/scripts/cluster_quality_check.sh -i $wd/maxbin \
    -o "$wd/maxbin" \
    -q $fastq \
    -x ".fasta"
  checkm_summary_to_tsv "$wd/maxbin/summary.txt"
  
  else "SKIP: $mes"
fi
# =================================================================================================
# VAMB


#vamb --outdir $wd/vamb --fasta $metagenome --bamfiles $bam -o C --minfasta 200000




cat "$wd/metabat2/summary.tsv" > checkm_all.tsv
tail -n+2 "$wd/maxbin/summary.txt" >> checkm_all.tsv





