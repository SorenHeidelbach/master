

######################################################################################################
# Multi to single
######################################################################################################

module load ont_fast5_api/1.4.7-foss-2018a-Python-3.6.4
# zymo WGA
multi_to_single_fast5 -i /scratch/users/sheide17/zymo_WGA_fast5_filt -s /scratch/users/sheide17/zymo_WGA_fast5_filt_single
# zymo NAT
multi_to_single_fast5 -i /scratch/users/sheide17/zymo_NAT_fast5_filt_100x -s /scratch/users/sheide17/zymo_NAT_fast5_filt_100x_single




######################################################################################################
# Tombo resquiggle
######################################################################################################

REF="/shared-nfs/SH/samples/zymo/assembly_100x.fasta"

# NAT
FASTQ_NAT="/shared-nfs/SH/samples/zymo/NAT/pooled_reads_filt_100x_trim.fastq"
FAST5_NAT="/scratch/users/sheide17/zymo_NAT_fast5_filt_100x_single"

# alignement (minimap2)
module load Minimap2/2.17-foss-2020b
module load SAMtools/1.11-foss-2020b
minimap2 -a -x map-ont $REF $FASTQ_NAT | samtools sort -T tmp -o $(dirname $REF)/output.sorted.bam
samtools index $(dirname $REF)/output.sorted.bam
# re-squggle (Tombo)
tombo resquiggle $FAST5_NAT $REF --processes 10 --corrected-group RawGenomeCorrected_001 --basecall-group Basecall_1D_000 --overwrite


# WGA
FASTQ_WGA="/shared-nfs/SH/samples/zymo/WGA/pooled_reads_filt_trim.fastq"
FAST5_WGA="/scratch/users/sheide17/zymo_WGA_fast5_filt_single"

# alignement (minimap2)
minimap2 -a -x map-ont $REF $FASTQ_WGA | samtools sort -T tmp -o $(dirname $REF)/output.sorted.bam
samtools index $(dirname $REF)/output.sorted.bam
# re-squggle (Tombo)
tombo resquiggle $FAST5_WGA $REF --processes 10 --corrected-group RawGenomeCorrected_002 --basecall-group Basecall_1D_000 --overwrite

