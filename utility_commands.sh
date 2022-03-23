

# Mapping for coverage
ref="zymoHMW_references.fasta"
fastq="PCR/basecalls_trim_filt.fastq"
file="PCR/trim_filt_ref"
module load Minimap2/2.17-foss-2020b
module load SAMtools/1.14-GCC-10.2.0
minimap2 -a $ref $fastq --secondary=no -x map-ont -t 50 |
    samtools sort -o $file.bam
module purge 

# Turn bam file into mappings summary
bam="/shared-nfs/SH/samples/zymoHMW/PCR/trim_filt_ref.bam"
file="trim_filt_ref"
module load SAMtools/1.14-GCC-10.2.0
# samtools view "$bam" > "$file.view.txt"
samtools depth --threads 20 -a "$bam" > "$file.depth"
module purge


