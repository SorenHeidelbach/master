
DIR="/user/student.aau.dk/sheide17/projects/current_difference/zymo"
cd $DIR
OUT="megalodon/pcr_test"
mkdir -p $OUT
FAST5="fast5_pcr_test"


REF="Zymo-Isolates-SPAdes-Illumina.fasta"
guppy_server_path="/user/student.aau.dk/sheide17/ont-guppy/bin/guppy_basecall_server"

megalodon \
  "/raid/student.sheide17/megalodon/zymo/fast5_pcr_test"  \
  --reference $REF \
  --output-directory $OUT \
  --outputs basecalls mappings signal_mappings \
  --guppy-config "dna_r9.4.1_450bps_hac.cfg" \
  --devices 0 \
  --guppy-server-path $guppy_server_path \
  --overwrite
