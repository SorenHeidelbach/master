


n = 8400
test <- signal_mappings_pcr_unnested %>% 
  filter(pos > n & pos < n+200 & direction == "rev") %>% 
  filter(read_id == read_id[1])
test[, dac_pos := as.numeric(pos) + ((1:.N)/.N), by = pos]
gg
test %>% 
  ggplot() +
  aes(x = dac_pos, y = dacs_norm) +
  #geom_vline(aes(xintercept = pos)) +
  geom_line(size = 0.5) +
  #geom_point(size = 1.3) +
  theme_void() +
  ylim(c(-2, 2)) +
  theme(axis.text = element_blank(), axis.title.y = element_blank()) +
  scale_y_continuous(n.breaks = 3) +
  labs(
    x = "Time"
  )
  
