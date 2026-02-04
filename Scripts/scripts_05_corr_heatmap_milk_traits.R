source(here("scripts", "00_setup.R"))

# build phyloseq for all matched samples
ps_species <- ps_base

make_metadata_corr_heatmap <- function(ps_obj, title_suffix) {
  
  if (ntaxa(ps_obj) > 50) {
    top_taxa <- names(sort(taxa_sums(ps_obj), decreasing = TRUE))[1:50]
    ps_obj <- prune_taxa(top_taxa, ps_obj)
  }
  
  otu_data <- as.matrix(otu_table(ps_obj))        # taxa x samples
  meta_sub <- data.frame(sample_data(ps_obj))
  
  target_cols <- c("Production","Fat","Energy","Lactose","Mineral","Protein")
  numeric_meta <- meta_sub[, target_cols, drop=FALSE]
  numeric_meta[] <- lapply(numeric_meta, as.numeric)
  
  corr_res <- psych::corr.test(t(otu_data), numeric_meta, method="spearman", adjust="none")
  cor_matrix <- corr_res$r
  p_matrix   <- corr_res$p
  
  star_matrix <- matrix("", nrow=nrow(p_matrix), ncol=ncol(p_matrix))
  star_matrix[p_matrix < 0.05]  <- "*"
  star_matrix[p_matrix < 0.01]  <- "**"
  star_matrix[p_matrix < 0.001] <- "***"
  
  my_color  <- colorRampPalette(c("#B2182B","white","#2166AC"))(100)
  my_breaks <- seq(-1, 1, length.out = 101)
  
  pheatmap(
    cor_matrix,
    main = paste(title_suffix, "Lactation Correlation"),
    color = my_color, breaks = my_breaks,
    display_numbers = star_matrix,
    fontsize_number = 14,
    cluster_rows = TRUE, cluster_cols = FALSE,
    fontsize_row = 10, fontsize_col = 12,
    angle_col = 45, border_color = "grey90"
  )
}

make_metadata_corr_heatmap(subset_samples(ps_species, Lactation == "Early"), "Early")
make_metadata_corr_heatmap(subset_samples(ps_species, Lactation == "Late"),  "Late")
