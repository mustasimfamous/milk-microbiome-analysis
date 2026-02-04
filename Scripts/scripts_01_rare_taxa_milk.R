source(here("scripts", "scripts_00_setup.R"))

# Milk samples only (if SampleType exists)
milk_ids <- meta_matched %>%
  filter(tolower(SampleType) == "milk") %>%
  pull(SampleID)

micro_milk <- micro_matched %>% select(Taxa, all_of(milk_ids))

# taxa x samples matrix
mat <- as.data.frame(micro_milk)
rownames(mat) <- mat$Taxa
mat$Taxa <- NULL

threshold <- 0.001
avg_abundance <- rowMeans(mat, na.rm = TRUE)
rare_taxa <- names(avg_abundance[avg_abundance < threshold])

rare_data <- mat[rare_taxa, , drop = FALSE]

top_20_rare <- rare_data %>%
  mutate(Mean = rowMeans(.)) %>%
  arrange(desc(Mean)) %>%
  head(20) %>%
  select(-Mean)

top_20_rare$Species <- rownames(top_20_rare)

plot_data <- top_20_rare %>%
  pivot_longer(cols = -Species, names_to = "SampleID", values_to = "Abundance")

p <- ggplot(plot_data, aes(x = SampleID, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  scale_fill_viridis_d(option = "plasma") +
  labs(
    title = "Top 20 Rare Taxa in Milk Samples",
    subtitle = paste("Threshold: Mean Abundance <", threshold),
    x = "Samples", y = "Relative Abundance", fill = "Rare Species"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8))

print(p)
ggsave(file.path(OUT_FIG, "01_rare_taxa_milk.png"), p, width = 14, height = 7, dpi = 600)
