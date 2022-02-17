#!/bin/bash
cd ~
# Main working directionary
MAINWD=/shared-nfs/SH
# Subfolder where everything will be done
SUBWD=zymo_mock_even
wd=$MAINWD/$SUBWD
mkdir $wd
cd $wd

## Location of scripts used
#SCRIPTS="/shared-nfs/RHK/supervision/bt6df20/testdata/SH/Scripts/"
# dir of data
READS="/shared-nfs/SH/zymo_mock_even/fastq/ERR3152364.fastq.gz"
# Filter size for filtlong
FILTER=0
FILTSIZE=500m
# Number of CPUs to use
THREADS=20
nRACON=0
CANUASSM=0 # Set assembly type 1 -> Canu, 0 -> Flye
GENOMESFASTANI=$MAINWD/Genomes # Folder of genomes to compare to bins with fastani
# Every tool action is saved to log and the time. > overwrites, >> apppends
date > log.txt

# Modules
MOD_BARRNAP=Barrnap/0.9-foss-2018a                      # https://github.com/tseemann/barrnap
MOD_CANU=canu/1.9-foss-2018a                            # https://readthedocs.org/projects/canu/downloads/pdf/latest/
MOD_CHECKM=CheckM/1.1.2-foss-2018a-Python-3.6.4         #
MOD_FASTANI=FastANI/1.2-foss-2018a                      # https://github.com/ParBLiSS/FastANI/blob/master/README.md
MOD_FILT=Filtlong/0.2.0-foss-2018a                      # https://github.com/rrwick/Filtlong
MOD_FLYE=Flye/2.8-1-foss-2018a-Python-3.6.4             # https://github.com/fenderglass/Flye
MOD_GTDBTK=GTDBTk/1.0.2-foss-2018a-Python-3.6.4         # https://github.com/Ecogenomics/GTDBTk
MOD_JAVA=Java/13.0.1
MOD_MAXBIN=MaxBin/2.2.7-foss-2018a-Perl-5.26.1          # https://sourceforge.net/projects/maxbin2/
MOD_MEDAKA=Medaka/1.2.3-foss-2020b                      # https://nanoporetech.github.io/medaka/installation.html#sequence-correction
MOD_METABAT=MetaBAT/2.12.1-foss-2018a                   # https://bitbucket.org/berkeleylab/metabat/src/master/
MOD_METAPHLAN=metaphlan2/2.9.21-foss-2018a              # https://github.com/brianmorganpalmer/metaphlan
MOD_MINIMAP=Minimap2/2.17-foss-2020b                    # https://github.com/lh3/minimap2
MOD_NANOPLOT=NanoPlot/1.24.0-foss-2018a                 # https://github.com/wdecoster/NanoPlot
MOD_PROKKA=prokka/1.14.0-foss-2018a-BioPerl-1.7.2       # https://github.com/tseemann/prokka
MOD_RACON=Racon/1.3.3-pikachu-foss-2018a                # https://github.com/lbcb-sci/racon
MOD_SAMTOOLS=SAMtools/1.10-foss-2018a                   # http://www.htslib.org/doc/samtools.html

######################################################################################################
#									          Filtering
######################################################################################################

if [$FILTER]; then \
	# Filtering, filtlong
	OUT=filteredreads.fastq
	TOOL="filtlong with $FILTSIZE bases"
	if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
		module load $MOD_FILT
		# Reduce number of bases in fastq to specified number
		filtlong --target_bases $FILTSIZE $READS > $OUT$FILTSIZE
		module purge
	if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
	READS="$wd/filteredreads.fastq"
´fi

# Quality check of filtered reads, Nanoplot
OUT=Nanoplot
TOOL=Nanoplot
if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
	module load $MOD_NANOPLOT
	NanoPlot -t $THREADS --verbose -o $OUT --fastq $READS
	module purge
if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi


######################################################################################################
#									          Assembly
######################################################################################################


if [$CANUASSM]; then \
	date +%X >> log.txt
	# Assembly, Canu
	OUT=Canu
	TOOL="Canu assembly"
	if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
	# Canu uses raw reads to make an assembly
	module load $MOD_CANU
	canu -p canu			    `# Set name of prefix to assembly`\
	   -d $OUT                  `# Set name of directory for assembly`\
	   -fast                    `# Speed up assemblt`\
	   genomeSize=5m            `# Set the size of the genome` \
	   -maxThreads=$THREADS\
	   -maxMemory=150\
	   -nanopore-raw $READS    `# Specify type of reads, here raw nanopore reads`
	if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
	date +%X >> log.txt
	ASSM="$wd/$OUT/canu.contigs.fasta"


else
	# Assembling, Flye
	OUT=Flye
	TOOL="flye with meta"
	if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
	module load $MOD_FLYE
	# Assembly with a guess of total contigs size of 30m and with metagenomic setting
	flye --nano-raw $READS --genome-size 30m --out-dir $OUT --meta --threads $THREADS
	module purge
	if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
	ASSM="$wd/$OUT/assembly.fasta"
	date +%X >> log.txt
fi

######################################################################################################
#									          Polishing
######################################################################################################


# Racon
n=1
while [[ $nRACON -ge $n ]]; do
	# Read alignment, minimap2
	# Used for racon polish. Alignment of ONT fastq reads to Flye assembly
	ALIGN=Minimap$n
	OUT=$ALIGN.sam
	TOOL="Minimap2"
	if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
	module load $MOD_MINIMAP
	module load $MOD_SAMTOOLS
	minimap2 -ax map-ont -t $THREADS $ASSM $READS > $OUT
	samtools view --threads $THREADS -Sb -F 0x104 $OUT |\
	samtools sort --threads $THREADS  -> $ALIGN.bam
	module purge
	if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
	date +%X >> log.txt
	# Racon
	TOOL="Racon"
	DIROUT=Racon$n
	if [ ! -s $DIROUT ]; then mkdir $DIROUT; fi
	OUT=$DIROUT/Racon.fasta
	if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
	module load $MOD_RACON
	racon -m 8 -x -6 -g -8 -w 500 -t $THREADS $READS $ALIGN.sam $ASSM > $OUT
	module purge
	if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
	ASSM="$wd/$OUT"
	$(( n++ ))
done
date +%X >> log.txt


# Medaka
OUT=Medaka
TOOL="Medaka"
if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
	module load $MOD_MEDAKA
	medaka_consensus -m r103_min_high_g345 -i $READS -d $ASSM -o $OUT -t $THREADS
	module purge
if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
ASSM="$wd/$OUT/consensus.fasta"
date +%X >> log.txt


######################################################################################################
#										 Binning
######################################################################################################


# Binning, Metabat
OUT=MetaBAT
BINS=$OUT
TOOL="Metabat"
if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
	module load $MOD_METABAT
	metabat --saveCls -i $ASSM -t $THREADS -o $OUT/bin
	module purge
if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
date +%X >> log.txt


# Check bin quality, checkM
OUT=CheckM
TOOL="checkM"
mkdir .temp
if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
module load $MOD_CHECKM
checkm lineage_wf $BINS $OUT -x fa -t $THREADS --tmpdir $wd/.temp --file summary.txt
checkm qa $OUT/lineage.ms $OUT > summary.txt
module purge
if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
date +%X >> log.txt


######################################################################################################
#										 Annotation
######################################################################################################


# Annotation, Prokka
OUT=Prokka
TOOL="Prokka"
if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
module load $MOD_PROKKA
module load $MOD_JAVA
for F in $BINS/*.fa; do
  N=$(basename $F .fa) ;
  prokka --locustag $N --outdir $OUT/$N --prefix $N --cpus $THREADS $F
done
module purge
if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
date +%X >> log.txt


# Generation of KO ID file to upload to KEGG reconstruct pathway
if [ -s $OUT/NameKO_ID.txt ]; then echo "Skipping KEPP ID generation" >> log.txt; else
rm $OUT/NameKO_ID.txt
rm $OUT/missedKO_ID.txt
for F in $OUT/*; do
  N=$(basename $F) ;
  echo "Working on $N"
  # Get KO IDs with names
  grep -o "UniProt.*" $OUT/$N/$N.gff | cut -d ";" -f1 | cut -d ":" -f2 > $OUT/$N/$N.uniprot.txt
  python3 $SCRIPTS/IDconvert.py $OUT/$N/$N.uniprot.txt "ACC+ID" "KO_ID" > $OUT/$N/$N.ID_KO.txt
  python3 $SCRIPTS/IDconvert.py $OUT/$N/$N.uniprot.txt "ACC+ID" "GENENAME" > $OUT/$N/$N.ID_name.txt
  Rscript $SCRIPTS/IDcombine.R $OUT/$N/$N.ID_name.txt $OUT/$N/$N.ID_KO.txt $OUT/namedKO_ID.txt $N >> $OUT/missedKO_ID.txt
done
fi


######################################################################################################
#										 Classification
######################################################################################################


# Classifying, GTDBTK
OUT=GTDBTk
TOOL="GTDBtk"
if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
module load $MOD_GTDBTK
gtdbtk classify_wf --genome_dir $BINS/ --out_dir $OUT --cpus $THREADS --extension fa
module purge
if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
date +%X >> log.txt


# rRNA, barrnap
OUT=Barrnap
TOOL="barrnap"
if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
mkdir $OUT
module load $MOD_BARRNAP
for F in $BINS/*.fa; do
  N=$(basename $F .fa) ;
  barrnap -o $OUT/$N.fa < $F > $OUT/$N.gff; done
module purge
if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi fi
date +%X >> log.txt









# Visual alignment with genomes of interest (GOI), fastANI
OUT=FastANI
TOOL="FastANI"
if [ -s $OUT ]; then echo "$OUT exist, skipping $TOOL" >> log.txt; else
module load $MOD_FASTANI
mkdir $OUT
for GOI in $GENOMESFASTANI/*; do
# Looping through reference genomes in genome folder
	GENOME=$(basename $GOI .fasta)
	mkdir $OUT/$GENOME
	for F in $BINS/*.fa; do
	  # Get basename of bin, ex. bin.1
	  N=$(basename $F .fa) ;
	  mkdir $OUT/$GENOME/$N
	  # Align with fastANI
	  fastANI -q $F -r $GOI --visualize --fragLen 2000 -o $OUT/$GENOME/$N/$N.out
	  # Rscript for visualization
	  Rscript $MAINWD/scripts/fastani_visualize.R $F $GOI $OUT/$GENOME/$N/$N.out.visual; done
done
if [ -s $OUT ]; then echo "Succesfully generated $OUT with $TOOL" >> log.txt; else echo "Failed to run $TOOL" >> log.txt; fi
module purge
fi
date +%X >> log.txt
echo "Done" >> log.txt
