source(here("scripts", "scripts_00_setup.R"))

keep_sample_name <- c(
  "142_AB","142_AB2","142_NI","142_NI2","180_AB","180_AB2","214_AB","214_AB2",
  "282_AB","282_AB2","354_AB","354_AB2","37_AB","37_AB2","490_NI","490_NI2",
  "503_NI","503_NI2","520_NI","520_NI2","563_NI","563_NI2","642_NI","642_NI2",
  "7358_NI","7358_NI2","7403_NI","7403_NI2","759_NI","759_NI2","797_NI","797_NI2",
  "989_NI","989_NI2"
)

data_filt <- meta_matched %>%
  filter(Sample_Name %in% keep_sample_name) %>%
  mutate(Stage = ifelse(grepl("(_AB|_NI)$", Sample_Name), "Early", "Late")) %>%
  mutate(Stage = factor(Stage, levels=c("Early","Late")))

trait_labels <- c(
  "Production"="Production (L/day)",
  "Energy"="Energy (kcal/L)",
  "Fat"="Fat (%)",
  "Protein"="Protein (%)",
  "Lactose"="Lactose (%)",
  "Mineral"="Minerals (%)"
)

long_data <- data_filt %>%
  select(Stage, all_of(names(trait_labels))) %>%
  pivot_longer(cols = all_of(names(trait_labels)),
               names_to="Trait", values_to="Value") %>%
  mutate(Trait_Full = recode(Trait, !!!trait_labels)) %>%
  mutate(Trait_Full = factor(Trait_Full, levels = trait_labels))

p <- ggplot(long_data, aes(Stage, Value, fill=Stage)) +
  geom_boxplot(width=0.6, alpha=0.9, color="black", linewidth=0.8) +
  facet_wrap(~Trait_Full, scales="free_y", ncol=3) +
  stat_compare_means(method="t.test", label="p.signif",
                     comparisons=list(c("Early","Late")), size=6) +
  theme_bw() +
  theme(legend.position="none",
        axis.text.x = element_text(face="bold", size=12))

print(p)
ggsave(file.path(OUT_FIG, "04_milk_traits_paired.png"), p, width = 14, height = 9, dpi = 600)

# summary table
summary_stats <- long_data %>%
  group_by(Trait, Stage) %>%
  summarise(Mean=mean(Value, na.rm=TRUE),
            SE=sd(Value, na.rm=TRUE)/sqrt(n()),
            .groups="drop")

p_vals <- long_data %>%
  group_by(Trait) %>%
  t_test(Value ~ Stage) %>%
  add_significance() %>%
  select(Trait, p, p.signif)

final_summary_table <- summary_stats %>%
  pivot_wider(names_from=Stage, values_from=c(Mean,SE)) %>%
  left_join(p_vals, by="Trait") %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

write.csv(final_summary_table, file.path(OUT_TAB, "04_milk_traits_summary.csv"), row.names = FALSE)
final_summary_table
