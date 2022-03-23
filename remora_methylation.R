pacman::p_load(
  "Rsamtools",
  "data.table",
  "ShortRead"
)

Rsamtools::sortBam(
  "/scratch/users/sheide17/samples/zymoHMW/PCR/basecalls_v2.bam",
  "/scratch/users/sheide17/samples/zymoHMW/PCR/basecalls_v2_sort"
)
bam_file <- Rsamtools::BamFile("/scratch/users/sheide17/samples/zymoHMW/PCR/basecalls_v2_sort.bam") 
bonito <- Rsamtools::scanBam(bam_file)



headers <- Rsamtools::scanBamHeader(bam_file)


asd <- fread("/scratch/users/sheide17/samples/zymoHMW/PCR/basecalls_v2_sort.view.txt")
