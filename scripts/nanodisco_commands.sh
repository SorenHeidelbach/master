nd_dir=/home/nanodisco/dataset
wd="$nd_dir$wd"

# Preprocessing
mes="Preprocess WGA"
if [[ ! -s $wd/preprocess/wga.sorted.bam ]]; then
  log "$mes"
  
  # WGA
  echo $nd_dir$fast5_nat_gzip
  nanodisco preprocess -p $threads \
    -f "$nd_dir$fast5_nat_gzip" \
    -o $wd/preprocess \
    -r $nd_dir$assembly \
    -s nat
  
  else log "SKIP: $mes"
fi

mes="Preprocess NAT "
if [[ ! -s $wd/preprocess/nat.sorted.bam ]]; then
  log "$mes"
  
  # NAT
  nanodisco preprocess -p $threads \
    -f $nd_dir$fast5_wga_gzip   \
    -o $wd/preprocess \
    -r $nd_dir$assembly \
    -s wga
  
  else log "SKIP: $mes"
fi

# Current difference
mes="Current differece"
if [[ ! -s $wd/difference ]]; then
  log "$mes"
  
  nanodisco difference -nj $threads -nc 10 -p 5 \
    -i $wd/preprocess \
    -o $wd/difference \
    -w wga \
    -n nat \
    -r $nd_dir$assembly
  
  else log "SKIP: $mes"
fi

nanodisco merge -d $wd/difference \
  -o $wd/difference_merge \
  -b sample

# Coverage
mes="Coverage WGA"
if [[ ! -s $wd/coverage/wga.cov ]]; then
  log "$mes"
  
  nanodisco coverage -b $wd/preprocess/zymo_wga.sorted.bam \
    -r $nd_dir$assembly \
    -o $wd/coverage
  
  else log "SKIP: $mes"
fi
mes="Coverage NAT"
if [[ ! -s $wd/coverage/nat.cov ]]; then
  log "$mes"
  
  nanodisco coverage -b $wd/preprocess/zymo_nat.sorted.bam \
    -r $nd_dir$assembly \
    -o $wd/coverage
  
  else log "SKIP: $mes"
fi

# Methylation profile
mes="Coverage NAT"
if [[ ! -s $wd/profile/methylation_profile_sample.RDS ]]; then
  log "$mes"
  
  nanodisco profile -p $threads \
    -r $nd_dir$assembly \
    -d $wd/difference_merge/zymo_difference.RDS \
    -w $wd/coverage/zymo_wga.cov \
    -n $wd/coverage/zymo_nat.cov \
    -b sample \
    -a "all" \
    -o $wd/profile \
    --min_contig_len 5000 \
    -c 5
  
  else log "SKIP: $mes"
fi
'
