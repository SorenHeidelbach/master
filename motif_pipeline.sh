conda activate py38_pyto

# Zymo
dif="/shared-nfs/SH/samples/zymo/nanodisco/difference_merge/zymo_difference.RDS"
out="/shared-nfs/SH/samples/zymo/modification_detection/mod_detection_7"
event_frame=10
threads=50
epochs=100
assembly="/shared-nfs/SH/samples/zymo/assembly_100x.fasta"

# MGM1
dif="/shared-nfs/SH/nanodisco/MGM1/dataset/metagenome_difference.RDS"
out="/shared-nfs/SH/samples/MGM1/mod_test1"
event_frame=10
threads=20
epochs=60
assembly="/shared-nfs/SH/nanodisco/MGM1/reference/metagenome.fasta"

mkdir -p $out
if [ ! -s "$out/event_features.tsv" ]; then
  mkdir $out
  motif_feature_preprocess.R \
    -d $dif \
    -o $out \
    --event_frame_size $event_frame \
    --event_u_val_threshold 0.00001
fi
if [ ! -s "$out/AE/$epochs.tsv" ]; then
  AE_modfind.py "$out/event_features.tsv" "$out/AE" \
    --event_frame_size $event_frame  \
    --epochs $epochs  \
    --batch_size 500 \
    --hidden_nodes 100 \
    --t $threads
fi

if [ ! -s "$out/DAAE/embedding_style_epoch_$epochs.tsv" ]; then
  DAAE_modfind.py "$out/event_features.tsv" "$out/DAAE" \
    --event_frame_size $event_frame  \
    --epochs $epochs  \
    --categorical_nodes 20 \
    --batch_size 500 \
    --hidden_nodes 100 \
    --t $threads
fi
if [ ! -s "$out/DAAE/UMAP_DBSCAN_consensus.pdf" ]; then
  motif_detection.R \
    -c  "$out/DAAE/embedding_categorical_epoch_$epochs.tsv" \
    -s  "$out/DAAE/embedding_style_epoch_$epochs.tsv" \
    -o  "$out" \
    -f  "$out/event_features.tsv" \
    -m  "$out/metainfo.tsv" \
    -a  $assembly \
    --run_umap \
    --train_loss "$out/DAAE/loss_training_epoch_$epochs.tsv" \
    --validation_loss "$out/DAAE/loss_validation_epoch_$epochs.tsv"
fi
cluster_quality_check.sh -i "$out/UMAP/bins" \
  -o "$out/UMAP/bins"
cluster_quality_check.sh -i "$out/DAAE/bins" \
  -o "$out/DAAE/bins"
# "/shared-nfs/SH/nanodisco/MGM1/reference/metagenome.fasta"



# Create bin membership  from bin fastas
echo -e "contig \t bin" > bin_membership.tsv
for i in *.fa; do
  awk -v file=${i%.fa} '/>/ { split($1, subfield, ">"); print subfield[2]"\t"file}' $i >> bin_membership.tsv
done
