# DEFINE GENE MODULES FOR ssGSEA
glycogen_modules <- list(
  
  Glycogen_Synthesis = c(
    "GYS1", "GYS2", "GBE1", "UGP2",
    "PGM1", "PPP1R3A", "PPP1R3B",
    "PPP1R3C", "PPP1R3D"
  ),
  
  Glycogen_Degradation = c(
    "PYGB", "PYGL", "PYGM",
    "AGL", "PHKA1", "PHKB", "PHKG1"
  ),
  
  Warburg_Pathway = c(
    "HK2", "LDHA", "PKM", "SLC2A1",
    "SLC2A3", "PFKL", "ENO1", "G6PD"
  ),
  
  Differentiation = c(
    "ESR1", "ESR2", "TP63", "KRT4",
    "KRT13", "KRT14", "IVL",
    "SPRR1A", "SPRR1B"
  )
)

# Show module sizes
cat("Gene modules defined:\n\n")
for (mod in names(glycogen_modules)) {
  cat(sprintf("  %-25s : %d genes\n",
              mod,
              length(glycogen_modules[[mod]])))
}



# COLLAPSE PROBES TO GENE LEVEL

# Add gene symbols to expression matrix rows
expr_with_symbols <- expr_matrix %>%
  as.data.frame() %>%
  mutate(PROBEID = rownames(.)) %>%
  left_join(probe_map, by = "PROBEID") %>%
  filter(!is.na(SYMBOL))

# Collapse: multiple probes per gene;use mean
expr_gene_level <- expr_with_symbols %>%
  dplyr::select(-PROBEID) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("SYMBOL")

cat("✓ Gene-level matrix built\n\n")
cat("  Original probes :", nrow(expr_matrix), "\n")
cat("  Unique genes now:", nrow(expr_gene_level), "\n")
cat("  Samples         :", ncol(expr_gene_level), "\n")




# RUN ssGSEA SCORING

library(GSVA)

# Version-safe ssGSEA execution
if (packageVersion("GSVA") >= "1.44") {
  
  # NEW syntax for GSVA >= 1.44
  ssgsea_param <- ssgseaParam(
    exprData  = as.matrix(expr_gene_level),
    geneSets  = glycogen_modules
  )
  
  ssgsea_scores <- gsva(ssgsea_param, verbose = FALSE)
  
} else {
  
  # OLD syntax for GSVA < 1.44
  ssgsea_scores <- gsva(
    expr          = as.matrix(expr_gene_level),
    gset.idx.list = glycogen_modules,
    method        = "ssgsea",
    kcdf          = "Gaussian",
    verbose       = FALSE
  )
}

# Check dimensions
cat("Score matrix dimensions:\n")
cat("  Modules :", nrow(ssgsea_scores), "\n")
cat("  Samples :", ncol(ssgsea_scores), "\n\n")

# first 5 sample scores
cat("Sample scores (first 5 samples):\n\n")
print(round(ssgsea_scores[, 1:5], 4))




# STEP 16a: CALCULATE GDI v1.0

# Extract individual module scores per sample
synthesis_score     <- ssgsea_scores["Glycogen_Synthesis", ]
degradation_score   <- ssgsea_scores["Glycogen_Degradation", ]
warburg_score       <- ssgsea_scores["Warburg_Pathway", ]
differentiation_score <- ssgsea_scores["Differentiation", ]

# Calculate GDI
# Higher numerator = more Warburg + more degradation
# Higher denominator = more synthesis + more differentiation
# Higher GDI = more depleted

GDI <- (warburg_score + degradation_score) /
  (synthesis_score + differentiation_score)

# Build results table
gdi_table <- data.frame(
  sample = names(GDI),
  GDI    = as.numeric(GDI),
  group  = sample_annotation$group
)

# Summary per group
cat("GDI mean per disease group:\n\n")

gdi_summary <- gdi_table %>%
  group_by(group) %>%
  summarise(
    n          = n(),
    mean_GDI   = round(mean(GDI),   4),
    median_GDI = round(median(GDI), 4),
    sd_GDI     = round(sd(GDI),     4),
    .groups    = "drop"
  )

print(gdi_summary)

# Quick trend check
normal_gdi <- gdi_summary$mean_GDI[gdi_summary$group == "Normal"]
cancer_gdi <- gdi_summary$mean_GDI[gdi_summary$group == "Cancer"]

cat("\n--- TREND CHECK ---\n")
if (cancer_gdi > normal_gdi) {
  cat("✓ GDI INCREASES from Normal to Cancer\n")
  cat("  Normal mean GDI :", normal_gdi, "\n")
  cat("  Cancer mean GDI :", cancer_gdi, "\n")
  cat("  Difference      :", round(cancer_gdi - normal_gdi, 4), "\n")
} else {
  cat("? Unexpected direction")
  cat("  Normal :", normal_gdi, "\n")
  cat("  Cancer :", cancer_gdi, "\n")
}




# STATISTICAL TESTING

# Are GDI values significantly different across groups?
kw_test <- kruskal.test(GDI ~ group, data = gdi_table)

cat("  Chi-squared:", round(kw_test$statistic, 3), "\n")
cat("  df         :", kw_test$parameter, "\n")
cat("  p-value    :", format(kw_test$p.value, scientific = TRUE), "\n")

if (kw_test$p.value < 0.05) {
  cat("  ✓ SIGNIFICANT — GDI differs across disease grades\n\n")
} else {
  cat(" Not significant")
}

# Does GDI increase monotonically with disease grade?
grade_numeric <- as.numeric(gdi_table$group)

spearman_test <- cor.test(
  gdi_table$GDI,
  grade_numeric,
  method = "spearman"
)

cat("Spearman Correlation (GDI vs disease grade):\n")
cat("  rho    :", round(spearman_test$estimate, 4), "\n")
cat("  p-value:", format(spearman_test$p.value, scientific = TRUE), "\n")

if (spearman_test$estimate > 0 & spearman_test$p.value < 0.05) {
  cat("  ✓ SIGNIFICANT positive correlation\n")
  cat("  GDI increases significantly with disease severity\n\n")
}

# Pairwise comparisons
# Which specific groups differ from Normal?
cat("Pairwise comparisons vs Normal (Wilcoxon test):\n\n")

normal_gdi_values <- gdi_table$GDI[gdi_table$group == "Normal"]
comparison_groups <- c("CIN1", "CIN2", "CIN3", "Cancer")

for (grp in comparison_groups) {
  grp_values <- gdi_table$GDI[gdi_table$group == grp]
  wt <- wilcox.test(grp_values, normal_gdi_values)
  sig <- ifelse(wt$p.value < 0.05, "✓ SIG", "— ns")
  cat(sprintf("  Normal vs %-8s : p = %s  %s\n",
              grp,
              format(wt$p.value, scientific = TRUE, digits = 3),
              sig))
}


# ROC CURVE ANALYSIS


# Install pROC if needed
if (!requireNamespace("pROC", quietly = TRUE)) {
  install.packages("pROC")
}
library(pROC)

cat("=== ROC CURVE ANALYSIS ===\n\n")

# Define binary outcomes
gdi_table <- gdi_table %>%
  mutate(
    CIN2plus = ifelse(group %in% c("CIN2","CIN3","Cancer"), 1, 0),
    CIN3plus = ifelse(group %in% c("CIN3","Cancer"), 1, 0)
  )

# ROC for CIN2+
roc_cin2 <- roc(gdi_table$CIN2plus, gdi_table$GDI, quiet = TRUE)
auc_cin2 <- round(auc(roc_cin2), 3)
ci_cin2  <- round(ci.auc(roc_cin2), 3)

# ROC for CIN3+
roc_cin3 <- roc(gdi_table$CIN3plus, gdi_table$GDI, quiet = TRUE)
auc_cin3 <- round(auc(roc_cin3), 3)
ci_cin3  <- round(ci.auc(roc_cin3), 3)

cat("AUC Results:\n\n")
cat(sprintf("  CIN2+ detection: AUC = %s (95%% CI: %s - %s)\n",
            auc_cin2, ci_cin2[1], ci_cin2[3]))
cat(sprintf("  CIN3+ detection: AUC = %s (95%% CI: %s - %s)\n",
            auc_cin3, ci_cin3[1], ci_cin3[3]))

# Interpret
cat("\nInterpretation:\n")
for (result in list(
  list(label="CIN2+", auc=auc_cin2),
  list(label="CIN3+", auc=auc_cin3)
)) {
  cat(sprintf("  GDI for %s → AUC = %s",
              result$label, result$auc))
  if (result$auc >= 0.80) {
    cat("  ✓ EXCELLENT discrimination\n")
  } else if (result$auc >= 0.70) {
    cat("  ✓ GOOD discrimination\n")
  } else if (result$auc >= 0.60) {
    cat("  ~ MODERATE discrimination\n")
  } else {
    cat(" POOR ")
  }
}




# FULL-RANGE GDI PLOT

# Make sure folder exists
dir.create("figures/GDI_plots", showWarnings = FALSE, recursive = TRUE)

# Full-range plot
gdi_plot_full <- ggplot(
  gdi_table,
  aes(x = group, y = GDI, fill = group)
) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.55,
    outlier.shape = 21,
    outlier.size = 2
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.7,
    alpha = 0.45,
    colour = "black"
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 5,
    colour = "white"
  ) +
  scale_fill_manual(values = group_colours) +
  labs(
    title = "Glycogen Depletion Index Across Cervical Disease Grades",
    subtitle = paste0(
      "Full-range plot showing all samples | ",
      "Kruskal-Wallis p = ",
      format(kw_test$p.value, scientific = TRUE, digits = 3),
      " | Spearman rho = ",
      round(spearman_test$estimate, 3)
    ),
    x = "Disease Grade",
    y = "Glycogen Depletion Index",
    caption = paste0(
      "GSE63514 | AUC CIN2+ = ", auc_cin2,
      " | AUC CIN3+ = ", auc_cin3,
      " | White diamond = group mean"
    )
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    plot.caption = element_text(size = 8, colour = "grey60"),
    axis.title = element_text(face = "bold"),
    legend.position = "none"
  )

# Save full-range plot
ggsave(
  filename = "figures/GDI_plots/GDI_v1_full_range.png",
  plot = gdi_plot_full,
  width = 10,
  height = 7,
  dpi = 300
)

print(gdi_plot_full)




# ZOOMED PUBLICATION GDI PLOT

# Use 95th percentile as upper zoom limit
# This prevents extreme cancer outliers from flattening the plot
y_upper <- quantile(gdi_table$GDI, 0.95, na.rm = TRUE)

cat("Zoomed y-axis upper limit:", round(y_upper, 3), "\n\n")

gdi_plot_zoom <- ggplot(
  gdi_table,
  aes(x = group, y = GDI, fill = group)
) +
  geom_boxplot(
    alpha = 0.75,
    width = 0.55,
    outlier.shape = 21,
    outlier.size = 2
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.8,
    alpha = 0.5,
    colour = "black"
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 23,
    size = 4.5,
    fill = "white",
    colour = "black"
  ) +
  scale_fill_manual(values = group_colours) +
  coord_cartesian(
    ylim = c(0, y_upper)
  ) +
  annotate(
    "text",
    x = 4.2,
    y = y_upper * 0.90,
    label = paste0(
      "CIN2+ AUC = ", auc_cin2,
      "\nCIN3+ AUC = ", auc_cin3
    ),
    size = 3.6,
    fontface = "bold",
    colour = "grey20"
  ) +
  labs(
    title = "Glycogen Depletion Index Increases Across Cervical Disease Grades",
    subtitle = paste0(
      "Higher GDI indicates greater glycogen depletion | ",
      "Kruskal-Wallis p = ",
      format(kw_test$p.value, scientific = TRUE, digits = 3),
      " | Spearman rho = ",
      round(spearman_test$estimate, 3)
    ),
    x = "Disease Grade",
    y = "Glycogen Depletion Index",
    caption = paste0(
      "GSE63514 | Zoomed to 95th percentile for readability; no data removed | ",
      "White diamond = median"
    )
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    plot.caption = element_text(size = 8, colour = "grey60"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    legend.position = "none"
  )

# Save zoomed plot
ggsave(
  filename = "figures/GDI_plots/GDI_v1_zoomed_publication.png",
  plot = gdi_plot_zoom,
  width = 10,
  height = 7,
  dpi = 300
)

print(gdi_plot_zoom)

cat("✓ Zoomed publication GDI plot saved:\n")
cat("  figures/GDI_plots/GDI_v1_zoomed_publication.png\n")



# ROC CURVE FIGURE

# Extract ROC curve coordinates
roc_cin2_df <- data.frame(
  specificity = roc_cin2$specificities,
  sensitivity = roc_cin2$sensitivities,
  group       = paste0("CIN2+  (AUC = ", auc_cin2, ")")
)

roc_cin3_df <- data.frame(
  specificity = roc_cin3$specificities,
  sensitivity = roc_cin3$sensitivities,
  group       = paste0("CIN3+  (AUC = ", auc_cin3, ")")
)

# Combine both curves
roc_combined <- rbind(roc_cin2_df, roc_cin3_df)

roc_colours <- setNames(
  c("#E74C3C", "#8E44AD"),
  c(
    paste0("CIN2+  (AUC = ", auc_cin2, ")"),
    paste0("CIN3+  (AUC = ", auc_cin3, ")")
  )
)

# Build ROC figure
roc_plot <- ggplot(
  roc_combined,
  aes(x = 1 - specificity,
      y = sensitivity,
      colour = group)
) +
  
  # Diagonal reference line (random chance)
  geom_abline(
    intercept = 0,
    slope     = 1,
    linetype  = "dashed",
    colour    = "grey60",
    size      = 0.5
  ) +
  
  # ROC curves
  geom_line(size = 1.2, alpha = 0.9) +
  
  # Colours
  scale_colour_manual(values = roc_colours) +
  
  # Axis formatting
  scale_x_continuous(
    limits = c(0, 1),
    labels = scales::percent_format()
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format()
  ) +
  
  # Labels
  labs(
    title    = "ROC Curves — Glycogen Depletion Index (GDI v1.0)",
    subtitle = paste0(
      "GDI as a classifier for clinically significant cervical dysplasia\n",
      "Dashed diagonal = random chance (AUC = 0.50)"
    ),
    x       = "1 - Specificity (False Positive Rate)",
    y       = "Sensitivity (True Positive Rate)",
    colour  = "Classification Target",
    caption = "GSE63514 | MolVIA Program — Publication 1"
  ) +
  
  theme_classic() +
  theme(
    plot.title      = element_text(face = "bold", size = 13),
    plot.subtitle   = element_text(size = 9, colour = "grey40"),
    plot.caption    = element_text(size = 8, colour = "grey60"),
    axis.title      = element_text(face = "bold"),
    legend.position = c(0.72, 0.18),
    legend.background = element_rect(
      fill    = "white",
      colour  = "grey80",
      size    = 0.4
    ),
    legend.title    = element_text(face = "bold", size = 9),
    legend.text     = element_text(size = 9)
  )

# Save
ggsave(
  filename = "figures/GDI_plots/GDI_v1_ROC_curves.png",
  plot     = roc_plot,
  width    = 8,
  height   = 7,
  dpi      = 300
)

print(roc_plot)

cat("✓ ROC curve figure saved:\n")
cat("  figures/GDI_plots/GDI_v1_ROC_curves.png\n")



# SAVE COMPLETE SESSION 3 WORKSPACE

# Save all objects
save(
  # Data objects
  expr_matrix,
  expr_gene_level,
  sample_annotation,
  probe_map,
  disease_order,
  group_colours,
  
  # Analysis objects
  design_matrix,
  contrast_matrix,
  fit1,
  fit2,
  
  # Results
  all_results,
  glycogen_results,
  glycogen_modules,
  
  # ssGSEA and GDI
  ssgsea_scores,
  gdi_table,
  gdi_summary,
  
  # Statistics
  kw_test,
  spearman_test,
  
  # ROC objects
  roc_cin2,
  roc_cin3,
  auc_cin2,
  auc_cin3,
  
  file = "data/processed/session03_complete.RData"
)
