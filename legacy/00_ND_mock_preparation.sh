# Download fast5 zipped
# Unzip fast5
# Extract fastq from fast5
# Random subsettnig with rasusa
# Extract read IDs from subsetted fastq
# Convert fast5 to multi fas5
# Subset fast5 with ont_fast5_api
export wd=/shared-nfs/SH/mock_ND
cd $wd

export links="https://sra-pub-src-2.s3.amazonaws.com/SRR10032543/MinION_HP_WGA_rep1.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032565/MinION_HP_WGA_rep2.tar.gz.1
https://sra-pub-src-1.s3.amazonaws.com/SRR10032544/MinION_HP_NAT.tar.gz.1

https://sra-pub-src-2.s3.amazonaws.com/SRR10032545/MinION_EC_WGA.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032546/MinION_EC_NAT.tar.gz.1

https://sra-pub-src-1.s3.amazonaws.com/SRR10032547/MinION_CP_WGA.tar.gz.1
https://sra-pub-src-1.s3.amazonaws.com/SRR10032548/MinION_CP_NAT.tar.gz.1

https://sra-pub-src-2.s3.amazonaws.com/SRR10032549/MinION_BF_WGA.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032555/MinION_BF_NAT.tar.gz.1

https://sra-pub-src-2.s3.amazonaws.com/SRR10032550/MinION_TP_WGA.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032551/MinION_TP_NAT.tar.gz.1

https://sra-pub-src-1.s3.amazonaws.com/SRR10032552/MinION_NO_WGA_2.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032553/MinION_NO_WGA_1.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032554/MinION_NO_NAT_2.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032556/MinION_NO_NAT_1.tar.gz.1

https://sra-pub-src-1.s3.amazonaws.com/SRR10032557/MinION_NG_WGA.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032558/MinION_NG_NAT.tar.gz.1

https://sra-pub-src-2.s3.amazonaws.com/SRR10032559/MinION_MH_WGA.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032560/MinION_MH_NAT.tar.gz.1

https://sra-pub-src-1.s3.amazonaws.com/SRR10032566/MinION_BA_WGA.tar.gz.1
https://sra-pub-src-2.s3.amazonaws.com/SRR10032567/MinION_BA_NAT.tar.gz.1"

export links="$links
MinION_HP_NAT_part1
MinION_HP_NAT_part2"


######################################################################################################
# Download
######################################################################################################
# Download file from each link
for val in $links; do
	name=$(echo $val | cut -d "/" -f 5)
	if [ ! -s /shared-nfs/SH/mock_ND/fast5_zipped/$name ]; then
	  echo $name
	  screen -d -m -S ${name} wget $val
	fi
done

######################################################################################################
# Unzip
######################################################################################################
# Unzip downloaded fast5
for val in $links; do
	name=$(echo $val | cut -d "/" -f 5 | cut -d "." -f 1)
	if [ ! -s /shared-nfs/SH/mock_ND/fast5/$name ]; then
	  echo $name
	  screen -d -m -S ${name} tar -xzf $name
	fi
done
# Checking status of unzipping
for val in $links; do
	name=$(echo $val | cut -d "/" -f 5 | cut -d "." -f 1)
	echo $name
	echo $(du -hs $name)
done

######################################################################################################
# Extracting Reads
######################################################################################################
# Adding manually unnested folder to "links"

screen -S "Extraction"
for val in $links; do
	export name=$(echo $val | cut -d "/" -f 5 | cut -d "." -f 1)
	if [ ! -s /shared-nfs/SH/mock_ND/fastq/${name}/${name}.fastq ]; then
  mkdir -p /shared-nfs/SH/mock_ND/fastq/$name
  module load nanopolish/0.11.3-foss-2018a-Python-3.6.4
  module purge
  module load nanopolish/0.11.3-foss-2018a-Python-3.6.4
  nanopolish extract ${wd}/fast5/${name}/* -o ${wd}/fastq/${name}/${name}.fastq
  module purge
	fi
done

######################################################################################################
# Random subsettnig with rasusa
######################################################################################################

conda activate py37_rasusa
for val in $links; do
	export name=$(echo $val | cut -d "/" -f 5 | cut -d "." -f 1)
	if [ ! -s /shared-nfs/SH/mock_ND/fastq_subset/${name}/${name}.fastq ]; then
  echo $name
  mkdir -p /shared-nfs/SH/mock_ND/fastq_subset/$name
  rasusa --coverage 100 --genome-size 3000000 --input ${wd}/fastq/${name}/${name}.fastq > \
  /shared-nfs/SH/mock_ND/fastq_subset/${name}/${name}.fastq
	fi
done

######################################################################################################
# Extract read IDs from subsetted fastq
######################################################################################################
for val in $links; do
  export name=$(echo $val | cut -d "/" -f 5 | cut -d "." -f 1)
  if [ ! -s ${wd}/readID/$name.txt ]; then
  echo $name
  mkdir -p ${wd}/readID
  awk '/^\@/ {print $1}' ${wd}/fastq_subset/${name}/${name}.fastq | \
  cut -d "_" -f 1 | \
  cut -d "@" -f 2 > ${wd}/readID/$name.txt
  fi
done

######################################################################################################
# Convert fast5 to multi fast5
######################################################################################################

screen -S "sinlge_to_multi"
module purge
module load ont_fast5_api/1.4.7-foss-2018a-Python-3.6.4
for val in $links; do
  export name=$(echo $val | cut -d "/" -f 5 | cut -d "." -f 1)
  if [ ! -s ${wd}/fast5_multi/$name ]; then
  echo $name
  single_to_multi_fast5 -i $wd/fast5/$name -s ${wd}/fast5_multi/$name --recursive -t 20
  fi
done

######################################################################################################
# Subset fast5 with ont_fast5_api
######################################################################################################

screen -S "fast5_subset"
module purge
module load ont_fast5_api/1.4.7-foss-2018a-Python-3.6.4
for val in $links; do
  export name=$(echo $val | cut -d "/" -f 5 | cut -d "." -f 1)
  if [ ! -s ${wd}/fast5_subset/$name ]; then
  echo $name
  fast5_subset -i ${wd}/fast5_multi/${name} -s ${wd}/fast5_subset/${name}/${name} -l ${wd}/readID/${name}.txt
  fi
done


######################################################################################################
# Pool subsetted fastq (only one WGA and one NAT from each organismn)
######################################################################################################
# KEEP: HP_WGA2, HP_NAT_part2, NO_NAT_1, NO_WGA_1
# DROP: 
drop="MinION_HP_WGA_rep1 
MinION_HP_NAT_part1 
MinION_NO_NAT_2 
MinION_NO_WGA_2"

rm $wd/fastq_pooled_WGA.fastq $wd/fastq_pooled_NAT.fastq
for path in $wd/fastq_subset/*; do
  name=$(basename $path)
  if ( ! ( echo $drop | grep -q $name ) ); then
    if [[ "$name" =~ WGA ]]; then
      echo "Writing to WGA"
      echo $name
      cat $path/$name.fastq >> $wd/fastq_pooled_WGA.fastq
    
    else
      echo "Writing to NAT"
      echo $name
      cat $path/$name.fastq >> $wd/fastq_pooled_NAT.fastq
    fi
  fi
done

######################################################################################################
# Pool subsetted fast5 (only one WGA and one NAT from each organismn)
######################################################################################################
# KEEP: HP_WGA2, HP_NAT_part2, NO_NAT_1, NO_WGA_1
# DROP: 
export drop="MinION_HP_WGA_rep1 
MinION_HP_NAT_part1 
MinION_NO_NAT_2 
MinION_NO_WGA_2"

screen -S pool_fast5_subset
mkdir $wd/fast5_pooled/WGA
mkdir $wd/fast5_pooled/NAT

for path in $wd/fast5_subset/*; do
  name=$(basename $path)
  if ( ! ( echo $drop | grep -q $name ) ); then
    if [[ "$name" =~ WGA ]]; then
      echo "Writing to WGA folder"
      echo $name
      cp $path/Min* $wd/fast5_pooled/WGA
    
    else
      echo "Writing to NAT folder"
      echo $name
      cp $path/Min* $wd/fast5_pooled/NAT
    fi
  fi
done



mkdir $wd/fast5_pooled/WGA_single
mkdir $wd/fast5_pooled/NAT_single

screen -S multi_to_single
module purge
module load ont_fast5_api/1.4.7-foss-2018a-Python-3.6.4
multi_to_single_fast5 -i $wd/fast5_pooled/WGA \
 -s $wd/fast5_pooled/WGA_single \
 -t 20
multi_to_single_fast5 -i $wd/fast5_pooled/NAT \
 -s $wd/fast5_pooled/NAT_single \
 -t 20


for path in $wd/fast5_pooled/WGA; do
  name=$(basename $path)
  
done





















