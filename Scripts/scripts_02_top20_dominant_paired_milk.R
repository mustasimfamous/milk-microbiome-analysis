source(here("scripts", "scripts_00_setup.R"))

keep_sample_name <- c(
  "142_AB","142_AB2","142_NI","142_NI2", "180_AB","180_AB2",
  "214_AB","214_AB2", "282_AB","282_AB2", "354_AB","354_AB2",
  "37_AB","37_AB2", "490_NI","490_NI2", "503_NI","503_NI2",
  "520_NI","520_NI2", "563_NI","563_NI2", "642_NI","642_NI2",
  "7358_NI","7358_NI2", "7403_NI","7403_NI2", "759_NI","759_NI2",
  "797_NI","797_NI2", "989_NI","989_NI2"
)

mapping <- meta_matched %>%
  filter(Sample_Name %in% keep_sample_name) %>%
  arrange(match(Sample_Name, keep_sample_name))

ids_to_keep <- intersect(mapping$SampleID, colnames(micro_matched)[-1])

micro_filtered <- micro_matched %>%
  select(Taxa, all_of(ids_to_keep))

mat <- as.data.frame(micro_filtered)
rownames(mat) <- mat$Taxa
mat$Taxa <- NULL
colnames(mat) <- mapping$Sample_Name[match(colnames(mat), mapping$SampleID)]

top_dominant <- mat %>%
  mutate(Mean = rowMeans(.)) %>%
  arrange(desc(Mean)) %>%
  head(20) %>%
  select(-Mean)

refine_names_unspecified <- function(tax_string) {
  p <- unlist(strsplit(tax_string, ";"))
  cl  <- if(length(p) >= 3) gsub("c__", "", p[3]) else "Unknown"
  ord <- if(length(p) >= 4) gsub("o__", "", p[4]) else "Unknown"
  fam <- if(length(p) >= 5) gsub("f__", "", p[5]) else ""
  gen <- if(length(p) >= 6) gsub("g__", "", p[6]) else ""
  spe <- if(length(p) >= 7) gsub("s__", "", p[7]) else ""
  
  if(!is.na(spe) && spe != "" && spe != "__" && nchar(spe) > 2) return(paste0("s__", spe))
  if(!is.na(gen) && gen != "" && gen != "__" && nchar(gen) > 2) return(paste0("g__", gen, " sp. (unspecified)"))
  if(!is.na(fam) && fam != "" && fam != "__" && nchar(fam) > 2) return(paste0("f__", fam, " gen. (unspecified)"))
  base_name <- if(ord != "Unknown" && ord != "") ord else cl
  paste0(base_name, " fam. (unspecified)")
}

top_dominant$TaxaName <- sapply(rownames(top_dominant), refine_names_unspecified)

plot_data <- top_dominant %>%
  pivot_longer(cols = -TaxaName, names_to = "Sample_Label", values_to = "Abundance")

plot_data$Sample_Label <- factor(plot_data$Sample_Label, levels = keep_sample_name)

p <- ggplot(plot_data, aes(x = Sample_Label, y = Abundance, fill = TaxaName)) +
  geom_bar(stat = "identity", width = 0.8) +
  theme_bw() +
  labs(title = "Top 20 Dominant Species Distribution",
       subtitle = "Unspecified naming for incomplete taxonomy",
       x = "Sample Name", y = "Relative Abundance", fill = "Microbes") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 10))

print(p)
ggsave(file.path(OUT_FIG, "02_top20_dominant_paired.png"), p, width = 16, height = 8, dpi = 600)
