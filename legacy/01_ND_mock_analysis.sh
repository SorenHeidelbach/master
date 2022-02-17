export wd=/shared-nfs/SH/mock_ND
cd $wd
export threads=20
export code=/shared-nfs/SH/code


######################################################################################################
# QC of sequencing runs
######################################################################################################
module load NanoPlot/1.38.0-foss-2020b
for path in $wd/fastq_subset/*; do
  name=$(basename $path)
  echo $name
  if [ ! -s ${wd}/NanoPlot/$name ]; then
  NanoPlot --fastq ${path}/${name}.fastq \
  -o ${wd}/NanoPlot/$name \
  -t $threads
  fi
done
module purge

######################################################################################################
# Assembly of mono cultures
######################################################################################################

# Medaka polish
screen -S medaka
module purge
module load Medaka/1.2.3-foss-2020b
for path in $wd/fastq_subset/*; do
  name=$(basename $path)
  if [ ! -s ${wd}/monoculture_assemblies/medaka/${name} ]; then
  medaka_consensus -m r941_min_high_g360 -i ${path}/${name}.fastq \
  -d ${wd}/monoculture_assemblies/racon3/${name}.fasta \
  -o ${wd}/monoculture_assemblies/medaka/${name} \
  -t $threads
  fi
done

######################################################################################################
# QC of mono culture assemblies
######################################################################################################

# Quast
mkdir -p $wd/monoculture/quast
module load QUAST/4.6.3-foss-2018a-Python-3.6.4
for path in ${wd}/monoculture/medaka/*; do
  name=$(basename $path)
  quast.py -o ${wd}/monoculture/quast/$name \
  -t $threads \
  ${wd}/monoculture/medaka/${name}/consensus.fasta
done

# CheckM
mkdir -p $wd/monoculture/assembly_temp
for path in ${wd}/monoculture/medaka/*; do
  name=$(basename $path)
  cat $path/consensus.fasta > $wd/monoculture/assembly_temp/$name.fasta
done

bash /shared-nfs/SH/code/bin_quality_check.sh -i $wd/monoculture/assembly_temp \
  -o $wd/monoculture/checkm \
  -x fasta \
  -t $threads
rm -r $wd/monoculture/assembly_temp

######################################################################################################
# Metagenomic assembly
######################################################################################################

# Assembly with flye
if [ ! -s ${wd}/metagenome/flye/assembly.fasta ]; then
module load Flye/2.9-GCC-10.2.0
flye --nano-raw $wd/fastq_pooled_NAT.fastq \
  --genome-size 30m \
  --out-dir $wd/metagenome/flye \
  --meta \
  --threads $threads
module purge
fi

# Polishing 
screen -S polishing
dos2unix $code/assembly_polish.sh
bash $code/fun_assembly_polish.sh \
  -r ${wd}/fastq_pooled_NAT.fastq \
  -a ${wd}/metagenome/flye/assembly.fasta \
  -o ${wd}/metagenome \
  -t $threads



######################################################################################################
# Nanodisco feature calculation
######################################################################################################

# Initiate screen
screen -S nanodisco

# Make sure Singularity is in PATH
echo 'export PATH=/shared-nfs/SH/software/singularity/bin:$PATH' >> ~/.bashrc && \
source ~/.bashrc

# Initiate Nanodisco container
singularity run --no-home -B /shared-nfs/SH/:/home/nanodisco/dataset -w /shared-nfs/SH/software/nanodisco/nd_env
nd_path=/home/nanodisco/dataset

# Preprocess (baseically nanopolish event alignment)
nanodisco preprocess -p 40 \
  -f $nd_path/mock_ND/fast5_pooled/NAT_single \
  -o $nd_path/mock_ND/metagenome/preprocess \
  -r $nd_path/mock_ND/metagenome/medaka/consensus.fasta \
  -s nd_mock_nat
nanodisco preprocess -p 40 \
  -f $nd_path/mock_ND/fast5_pooled/WGA_single \
  -o $nd_path/mock_ND/metagenome/preprocess \
  -r $nd_path/mock_ND/metagenome/medaka/consensus.fasta \
  -s nd_mock_wga

# Current differences

nanodisco difference -nj 20 -nc 4 -p 10 -f 1 -l 6973  \
  -i $nd_path/mock_ND/metagenome/preprocess \
  -o $nd_path/mock_ND/metagenome/difference \
  -w nd_mock_wga \
  -n nd_mock_nat \
  -r $nd_path/mock_ND/metagenome/medaka/consensus.fasta


# Calculate methylation profile







