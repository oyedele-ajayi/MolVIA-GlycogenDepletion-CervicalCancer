# PROBE TO GENE MAPPING

library(hgu133plus2.db)
library(AnnotationDbi)

# Get mapping table
probe_map <- select(
  hgu133plus2.db,
  keys    = rownames(expr_matrix),
  columns = c("PROBEID", "SYMBOL"),
  keytype = "PROBEID"
)

# Remove probes with no gene symbol
probe_map <- probe_map %>% filter(!is.na(SYMBOL))

cat("  Probes with gene symbol:", nrow(probe_map), "\n")
cat("  Unique genes represented:", length(unique(probe_map$SYMBOL)), "\n\n")

# Save for later use
write.csv(probe_map, "data/processed/probe_to_gene.csv", row.names = FALSE)


# SEARCH FOR KEY GLYCOGEN GENES

key_genes <- c(
  "GYS1", "GYS2", "UGP2", "GBE1", "PPP1R3C", "PPP1R3D",
  "PYGB", "PYGL", "GSK3B", "AKT1",
  "HK2", "LDHA", "PKM", "ESR1", "KRT4", "KRT13"
)

# SEARCHING FOR KEY GENES
gene_check <- probe_map %>%
  filter(SYMBOL %in% key_genes) %>%
  group_by(SYMBOL) %>%
  summarise(probe_count = n(), .groups = "drop") %>%
  arrange(SYMBOL)
print(gene_check)

found   <- nrow(gene_check)
missing <- setdiff(key_genes, gene_check$SYMBOL)

cat("\n✓ Found:", found, "of", length(key_genes), "genes\n")
if (length(missing) > 0) {
  cat("⚠ Missing:", paste(missing, collapse = ", "), "\n")
} else {
  cat("✓ All key genes present!\n")
}


# GYS1 AND HK2 EXPRESSION PLOTS

# Helper function to plot a single gene
plot_gene_expression <- function(gene_name, expr_mat, samp_annot, probe_mapping) {
  
  # Find probe for this gene (take first if multiple)
  probe_id <- probe_mapping %>%
    filter(SYMBOL == gene_name) %>%
    pull(PROBEID) %>%
    .[1]
  
  if (is.na(probe_id)) {
    return(NULL)
  }
  
  # Build data frame
  gene_df <- data.frame(
    expression = as.numeric(expr_mat[probe_id, ]),
    group      = samp_annot$group
  )
  
  # Calculate means
  means <- gene_df %>%
    group_by(group) %>%
    summarise(mean_expr = round(mean(expression), 3), .groups = "drop")
  
  cat("=== ", gene_name, " MEAN EXPRESSION ===\n")
  print(means)
  cat("\n")
  
  # Create plot
  p <- ggplot(gene_df, aes(x = group, y = expression, fill = group)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = 21, outlier.size = 2) +
    geom_jitter(width = 0.12, size = 1.5, alpha = 0.4, colour = "black") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, colour = "white") +
    scale_fill_manual(values = group_colours) +
    labs(
      title = paste0(gene_name, " Expression Across Cervical Disease Grades"),
      subtitle = "White diamond = group mean",
      x = "Disease Grade", y = "Expression (log2)"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      legend.position = "none"
    )
  
  return(p)
}

# Generate plots
gys1_plot <- plot_gene_expression("GYS1", expr_matrix, sample_annotation, probe_map)
hk2_plot  <- plot_gene_expression("HK2",  expr_matrix, sample_annotation, probe_map)

# Save plots
if (!is.null(gys1_plot)) {
  ggsave("figures/quality_control/QC_04_GYS1.png", gys1_plot, width = 8, height = 6, dpi = 300)
  print(gys1_plot)
}

if (!is.null(hk2_plot)) {
  ggsave("figures/quality_control/QC_05_HK2.png", hk2_plot, width = 8, height = 6, dpi = 300)
  print(hk2_plot)
}



# GYS1 AND HK2 EXPRESSION PLOTS Combine plot

# Load patchwork for combining plots
library(patchwork)

# Helper function to plot a single gene
plot_gene_expression <- function(gene_name, expr_mat, samp_annot, probe_mapping) {
  
  # Find probe for this gene (take first if multiple)
  probe_id <- probe_mapping %>%
    filter(SYMBOL == gene_name) %>%
    pull(PROBEID) %>%
    .[1]
  
  if (is.na(probe_id)) {
    return(NULL)
  }
  
  # Build data frame
  gene_df <- data.frame(
    expression = as.numeric(expr_mat[probe_id, ]),
    group      = samp_annot$group
  )
  
  # Calculate means
  means <- gene_df %>%
    group_by(group) %>%
    summarise(mean_expr = round(mean(expression), 3), .groups = "drop")
  
  cat("=== ", gene_name, " MEAN EXPRESSION ===\n")
  print(means)
  cat("\n")
  
  # Create plot
  p <- ggplot(gene_df, aes(x = group, y = expression, fill = group)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = 21, outlier.size = 2) +
    geom_jitter(width = 0.12, size = 1.5, alpha = 0.4, colour = "black") +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, colour = "white") +
    scale_fill_manual(values = group_colours) +
    labs(
      title    = paste0(gene_name, " Expression Across Cervical Disease Grades"),
      subtitle = "White diamond = group mean",
      x        = "Disease Grade",
      y        = "Expression (log2)"
    ) +
    theme_classic() +
    theme(
      plot.title    = element_text(face = "bold", size = 12),
      legend.position = "none"
    )
  
  return(p)
}

# Generate individual plots
gys1_plot <- plot_gene_expression("GYS1", expr_matrix, sample_annotation, probe_map)
hk2_plot  <- plot_gene_expression("HK2",  expr_matrix, sample_annotation, probe_map)

# Save individual plots
if (!is.null(gys1_plot)) {
  ggsave("figures/quality_control/QC_04_GYS1.png", gys1_plot, width = 8, height = 6, dpi = 300)
}

if (!is.null(hk2_plot)) {
  ggsave("figures/quality_control/QC_05_HK2.png", hk2_plot, width = 8, height = 6, dpi = 300)
}

# Combine plots side by side using patchwork
if (!is.null(gys1_plot) && !is.null(hk2_plot)) {
  
  combined_plot <- gys1_plot + hk2_plot +
    plot_layout(ncol = 2) +
    plot_annotation(
      # Overall title and subtitle for the combined figure
      title    = "Glycogen Metabolism Gene Expression Across Cervical Disease Grades",
      subtitle = "Left: GYS1 (Glycogen Synthase 1)  |  Right: HK2 (Hexokinase 2)",
      theme = theme(
        plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5)
      )
    )
  
  # Save combined figure
  ggsave(
    "figures/quality_control/QC_04_05_GYS1_HK2_combined.png",
    combined_plot,
    width  = 16,
    height = 6,
    dpi    = 300
  )
  
  # Display combined figure
  print(combined_plot)
  
  cat("✓ Combined GYS1 + HK2 side-by-side plot saved\n")
}


# SAVE EVERYTHING FOR NEXT SESSION

save(
  expr_matrix,
  sample_annotation,
  sample_info,
  probe_map,
  pca_result,
  disease_order,
  group_colours,
  file = "data/processed/session02_complete.RData"
)
