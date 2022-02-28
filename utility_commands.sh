

# Turn bam file into mappings summary
file=mappings
module load SAMtools/1.14-GCC-10.2.0
samtools sort $file.bam |
    samtools view > $file.sort.view.txt