#!/bin/bash

######################################################################################################
# Initialisation
######################################################################################################

# Help function for inputs
helpFunction()
{
 echo ""
 echo "Usage: $0 -r Reads -a Assembly -o Output_dir-t [Threads]"
 echo -r " read FASTQ file(s)"
 echo -a " assembly that need polishing"
 echo -o " output directionary of all racon and medaka runs"
 echo -t " number of threads to use"
 exit 1 # Exit script after printing help
}
threads=10

# Reads input options
while getopts "r:a:o:t:" opt; do
  case "$opt" in
  h) helpFunction ;; # Print helpFunction in case parameter is non-existent
  r) read_fastq="$OPTARG" ;;
  a) assembly="$OPTARG" ;;
  o) o_dir="$OPTARG" ;;
  t) threads="$OPTARG" ;;
  esac
done
echo test
# Check for essential input options
shift "$(( OPTIND - 1 ))"
if [ -z "$read_fastq" ] || [ -z "$assembly" ] || [ -z "$o_dir" ]; then
  helpFunction
fi

######################################################################################################
# Function
######################################################################################################
# This script takes in an assembly and accompaniyng fastq read files
# It runs 3x Racon and 1x medaka
echo down_to_actual_script
mkdir -p $o_dir
module load Racon/1.3.3-pikachu-foss-2018a
module load Minimap2/2.17-foss-2020b

# Racon 1
if [ ! -s $o_dir/racon1/assembly.fasta ]; then
echo "Racon 1" >> $o_dir/log.txt
mkdir -p $o_dir/racon1
minimap2 -x map-ont -t $threads \
  $assembly \
  $read_fastq > \
  $o_dir/racon1/alignment.paf
  
racon -t $threads \
  $read_fastq \
  $o_dir/racon1/alignment.paf \
  $assembly > \
  $o_dir/racon1/assembly.fasta 
else 
echo "$o_dir/racon1/assembly.fasta was identified and therefore racon1 skipped" >> $o_dir/log.txt
fi

# Racon 2
if [ ! -s $o_dir/racon2/assembly.fasta ]; then
echo "Racon 2" >> $o_dir/log.txt
mkdir -p $o_dir/racon2
minimap2 -x map-ont -t $threads \
  $o_dir/racon1/assembly.fasta  \
  $read_fastq > \
  $o_dir/racon2/alignment.paf
  
racon -t $threads \
  $read_fastq \
  $o_dir/racon2/alignment.paf \
  $o_dir/racon1/assembly.fasta > \
  $o_dir/racon2/assembly.fasta 
else 
echo "$o_dir/racon2/assembly.fasta was identified and therefore racon2 skipped" >> $o_dir/log.txt
fi

# Racon 3
if [ ! -s $o_dir/racon3/assembly.fasta ]; then
echo "Racon 3" >> $o_dir/log.txt
mkdir -p $o_dir/racon3
minimap2 -x map-ont -t $threads \
  $o_dir/racon2/assembly.fasta  \
  $read_fastq > \
  $o_dir/racon3/alignment.paf
  
racon -t $threads \
  $read_fastq \
  $o_dir/racon3/alignment.paf \
  $o_dir/racon2/assembly.fasta > \
  $o_dir/racon3/assembly.fasta 
else 
echo "$o_dir/racon3/assembly.fasta was identified and therefore racon3 skipped" >> $o_dir/log.txt
fi
module purge


module purge
module load Medaka/1.2.3-foss-2020b
if [ ! -s $o_dir/medaka/ ]; then
echo "Medaka" >> $o_dir/log.txt
medaka_consensus -m r941_min_high_g360 \
  -i $read_fastq \
  -d $o_dir/racon3/assembly.fasta \
  -o $o_dir/medaka/ \
  -t $threads
else 
echo "$o_dir/medaka/assembly.fasta was identified and therefore medaka skipped" >> $o_dir/log.txt
fi


