source(here("scripts", "scripts_00_setup.R"))

keep_sample_name <- c(
  "142_AB","142_AB2","142_NI","142_NI2","180_AB","180_AB2","214_AB","214_AB2",
  "282_AB","282_AB2","354_AB","354_AB2","37_AB","37_AB2","490_NI","490_NI2",
  "503_NI","503_NI2","520_NI","520_NI2","563_NI","563_NI2","642_NI","642_NI2",
  "7358_NI","7358_NI2","7403_NI","7403_NI2","759_NI","759_NI2","797_NI","797_NI2",
  "989_NI","989_NI2"
)

ps <- subset_samples(ps_base,
                     Sample_Name %in% keep_sample_name &
                       tolower(Lactation) %in% c("early","late"))
ps <- prune_samples(sample_names(ps), ps)

# relative abundance (%)
ps_rel <- transform_sample_counts(ps, function(x) x/sum(x) * 100)

abund <- t(as(otu_table(ps_rel), "matrix"))  # samples x taxa
top_taxa <- names(sort(colSums(abund), decreasing = TRUE))[1:min(20, ncol(abund))]

abund_plot <- cbind(
  abund[, top_taxa, drop=FALSE],
  Other = rowSums(abund[, setdiff(colnames(abund), top_taxa), drop=FALSE])
)

plot_df <- as.data.frame(abund_plot) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta_matched %>% select(SampleID, Lactation), by="SampleID") %>%
  mutate(Lactation = ifelse(tolower(Lactation)=="early","Early","Late")) %>%
  pivot_longer(cols = -c(SampleID, Lactation), names_to="Taxa", values_to="RelAbund")

sum_df <- plot_df %>%
  group_by(Lactation, Taxa) %>%
  summarise(RelAbund = mean(RelAbund), .groups="drop") %>%
  mutate(Lactation = factor(Lactation, levels=c("Early","Late")))

p <- ggplot(sum_df, aes(Lactation, RelAbund, fill = Taxa)) +
  geom_col(width = 0.7, color="white", linewidth=0.25) +
  theme_bw(base_size = 14) +
  labs(x=NULL, y="Mean relative abundance (%)", fill="Top taxa") +
  theme(panel.grid.major.x = element_blank())

print(p)
ggsave(file.path(OUT_FIG, "03_dominant_early_vs_late.png"), p, width = 8, height = 7, dpi = 600)
