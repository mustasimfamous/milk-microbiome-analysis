# scripts_00_setup.R
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(phyloseq)
  library(vegan)
  library(rstatix)
  library(ggpubr)
  library(pheatmap)
  library(psych)
  library(patchwork)
  library(here)
})

# paths
DATA_DIR <- here("data")
OUT_FIG  <- here("outputs", "figures")
OUT_TAB  <- here("outputs", "tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)

# read data ONCE
micro_raw <- read.delim(
  file.path(DATA_DIR, "species_rel_no30.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE
)
colnames(micro_raw)[1] <- "Taxa"

meta <- read.csv(
  file.path(DATA_DIR, "metadata_clean.csv"),
  stringsAsFactors = FALSE
) %>%
  mutate(
    SampleID   = str_trim(as.character(SampleID)),
    SampleType = str_trim(as.character(SampleType)),
    Lactation  = str_trim(as.character(Lactation)),
    Sample_Name = if ("Sample_Name" %in% names(.)) str_trim(as.character(Sample_Name)) else Sample_Name
  )

# match samples ONCE
otu_samples <- colnames(micro_raw)[-1]
common_samples <- intersect(otu_samples, meta$SampleID)
if (length(common_samples) == 0) stop("No matching SampleID between OTU table and metadata")

micro_matched <- micro_raw %>% select(Taxa, all_of(common_samples))
meta_matched  <- meta %>% filter(SampleID %in% common_samples) %>% slice(match(common_samples, SampleID))

# build base phyloseq ONCE
otu_mat <- as.matrix(micro_matched[,-1])
rownames(otu_mat) <- micro_matched$Taxa
mode(otu_mat) <- "numeric"

ps_base <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  sample_data({ df <- meta_matched; rownames(df) <- df$SampleID; df })
)
