source(here("scripts", "scripts_00_setup.R"))

# milk only early/late
ps_milk <- subset_samples(ps_base,
                          tolower(SampleType) == "milk" &
                            Lactation %in% c("Early","Late"))
ps_milk <- prune_samples(sample_names(ps_milk), ps_milk)

print(table(sample_data(ps_milk)$Lactation))

# alpha from relative abundance (Observed/Shannon/Simpson)
otu_m <- t(as(otu_table(ps_milk), "matrix"))
otu_pa <- otu_m; otu_pa[otu_pa > 0] <- 1

alpha_df <- data.frame(
  SampleID = rownames(otu_m),
  Observed = specnumber(otu_pa),
  Shannon  = diversity(otu_m, index="shannon"),
  Simpson  = diversity(otu_m, index="simpson"),
  stringsAsFactors = FALSE
) %>%
  left_join(data.frame(sample_data(ps_milk)) %>% select(SampleID, Lactation),
            by="SampleID") %>%
  pivot_longer(cols=c(Observed,Shannon,Simpson), names_to="Measure", values_to="Value") %>%
  mutate(Lactation = factor(Lactation, levels=c("Early","Late")))

lact_cols <- c("Early"="#1F78B4", "Late"="#E7298A")

alpha_plot_sig <- function(measure, ylab) {
  dfm <- alpha_df %>% filter(Measure == measure)
  y_max <- max(dfm$Value, na.rm=TRUE) * 1.1
  
  stat <- dfm %>%
    wilcox_test(Value ~ Lactation) %>%
    mutate(
      p.adj = p.adjust(p, method="BH"),
      p.signif = case_when(
        p.adj <= 0.001 ~ "***",
        p.adj <= 0.01  ~ "**",
        p.adj <= 0.05  ~ "*",
        TRUE ~ "ns"
      ),
      group1="Early", group2="Late",
      y.position = y_max
    ) %>% filter(p.adj <= 0.05)
  
  p <- ggplot(dfm, aes(Lactation, Value, fill=Lactation)) +
    geom_boxplot(width=0.65, alpha=0.9, color="black") +
    scale_fill_manual(values=lact_cols) +
    theme_classic(base_size=16) +
    labs(title=measure, x="Lactation", y=ylab)
  
  if (nrow(stat) > 0) {
    p <- p + stat_pvalue_manual(stat, label="p.signif",
                                xmin="group1", xmax="group2",
                                y.position="y.position",
                                inherit.aes=FALSE)
  }
  p
}

p_obs <- alpha_plot_sig("Observed","Observed richness")
p_sha <- alpha_plot_sig("Shannon","Shannon index")
p_sim <- alpha_plot_sig("Simpson","Simpson index")

pcoa_plot <- function(ps, dist_method, title_text) {
  ord <- ordinate(ps, "PCoA", dist_method)
  df  <- plot_ordination(ps, ord, justDF=TRUE) %>%
    mutate(Lactation = factor(Lactation, levels=c("Early","Late")))
  df_cloud <- df %>% group_by(Lactation) %>% filter(n() >= 3) %>% ungroup()
  
  p1 <- round(ord$values$Relative_eig[1]*100, 1)
  p2 <- round(ord$values$Relative_eig[2]*100, 1)
  
  ggplot(df, aes(Axis.1, Axis.2)) +
    stat_ellipse(data=df_cloud, aes(group=Lactation, fill=Lactation),
                 geom="polygon", alpha=0.20, color=NA) +
    stat_ellipse(data=df_cloud, aes(group=Lactation, color=Lactation),
                 linewidth=0.6) +
    geom_point(aes(color=Lactation), shape=17, size=4.5, alpha=0.95) +
    scale_color_manual(values=lact_cols) +
    scale_fill_manual(values=lact_cols) +
    theme_classic(base_size=16) +
    labs(title=title_text,
         x=paste0("PCoA1 (",p1,"%)"),
         y=paste0("PCoA2 (",p2,"%)"),
         color="Lactation") +
    guides(fill="none")
}

p_bray <- pcoa_plot(ps_milk, "bray", "Bray–Curtis PCoA (Milk)")
p_jacc <- pcoa_plot(ps_milk, "jaccard", "Jaccard PCoA (Milk)")

final_6panel <- (p_obs | p_sha) /
  (p_sim | p_bray) /
  (p_jacc | plot_spacer())

final_6panel
ggsave(file.path(OUT_FIG, "06_milk_alpha_beta_6panel.png"), final_6panel, width=14, height=15, dpi=600)
