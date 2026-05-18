# LOAD SESSION 2 DATA

# Clear environment for clean start
rm(list = ls())

# Load all saved objects from previous session
load("data/processed/session02_complete.RData")

# Load libraries required
library(tidyverse)
library(limma)
library(ggplot2)
library(ggrepel)


# Confirm data loaded correctly
cat("Expression matrix:", nrow(expr_matrix), "probes x",
    ncol(expr_matrix), "samples\n")
cat("Sample annotation:", nrow(sample_annotation), "samples\n\n")
cat("Group counts:\n")
print(table(sample_annotation$group))



# VERIFY SAMPLE ORDER IS CONSISTENT


# Check if column names of matrix match annotation titles
match_check <- all(colnames(expr_matrix) == sample_annotation$title)

if (match_check) {
  cat("✓ Sample order is CONSISTENT\n")
  cat("  Expression matrix columns match annotation rows exactly\n")
} else {
  cat("MISMATCH DETECTED — fixing now\n")
  
  # Reorder annotation to match matrix columns
  sample_annotation <- sample_annotation[
    match(colnames(expr_matrix), sample_annotation$title), ]
  
  # Verify fix worked
  if (all(colnames(expr_matrix) == sample_annotation$title)) {
    cat("Mismatch fixed successfully\n")
  } else {
    cat("Could not be fixed\n")
  }
}

cat("\nFirst 5 matrix columns:", colnames(expr_matrix)[1:5], "\n")
cat("First 5 annotation rows:", sample_annotation$title[1:5], "\n")



# BUILD LIMMA DESIGN MATRIX


# Build design matrix from group factor
design_matrix <- model.matrix(~ 0 + group, data = sample_annotation)

# Clean up column names — remove "group" prefix
colnames(design_matrix) <- gsub("group", "", colnames(design_matrix))

print(head(design_matrix, 10))

cat("\nDesign matrix dimensions:\n")
cat("  Rows (samples):", nrow(design_matrix), "\n")
cat("  Cols (groups) :", ncol(design_matrix), "\n")

cat("\nColumn sums (= samples per group):\n")
print(colSums(design_matrix))



# DEFINE STATISTICAL CONTRASTS

contrast_matrix <- makeContrasts(
  
  # Primary comparisons: each grade vs Normal
  CIN1_vs_Normal   = CIN1   - Normal,
  CIN2_vs_Normal   = CIN2   - Normal,
  CIN3_vs_Normal   = CIN3   - Normal,
  Cancer_vs_Normal = Cancer - Normal,
  
  # Precancer progression comparison
  CIN3_vs_CIN1     = CIN3   - CIN1,
  
  levels = design_matrix
)

cat("=== CONTRAST MATRIX ===\n\n")
print(contrast_matrix)

cat("\n✓ We will test", ncol(contrast_matrix), "comparisons:\n")
cat("  1. CIN1 vs Normal\n")
cat("  2. CIN2 vs Normal\n")
cat("  3. CIN3 vs Normal\n")
cat("  4. Cancer vs Normal\n")
cat("  5. CIN3 vs CIN1 (precancer progression)\n")



# FIT LINEAR MODEL TO ALL PROBES

cat("Fitting linear model to", nrow(expr_matrix), "probes...\n")


# Step 1: Fit model to each probe
fit1 <- lmFit(
  object = expr_matrix,
  design = design_matrix
)

cat("  Probes modelled:", nrow(fit1$coefficients), "\n")
cat("  Groups modelled:", ncol(fit1$coefficients), "\n")




# APPLY CONTRASTS AND EMPIRICAL BAYES

# Apply our contrasts
fit2 <- contrasts.fit(fit1, contrast_matrix)

# Apply empirical Bayes moderation
fit2 <- eBayes(fit2)


# Quick summary of significant genes per contrast
cat("Significant genes (FDR < 0.05, |logFC| > 1):\n\n")

for (contrast_name in colnames(contrast_matrix)) {
  
  results <- topTable(
    fit2,
    coef       = contrast_name,
    number     = Inf,
    adjust.method = "BH",
    sort.by    = "P"
  )
  
  sig_count <- sum(results$adj.P.Val < 0.05 & abs(results$logFC) > 1)
  up_count  <- sum(results$adj.P.Val < 0.05 & results$logFC > 1)
  down_count <- sum(results$adj.P.Val < 0.05 & results$logFC < -1)
  
  cat(sprintf("  %-20s: %d significant (%d up, %d down)\n",
              contrast_name, sig_count, up_count, down_count))
}




# EXTRACT AND SAVE FULL RESULTS TABLES

# Store all results in a list
all_results <- list()

for (contrast_name in colnames(contrast_matrix)) {
  
  # Extract full results table
  res <- topTable(
    fit2,
    coef          = contrast_name,
    number        = Inf,
    adjust.method = "BH",
    sort.by       = "P"
  )
  
  # Add probe ID as column
  res$PROBEID <- rownames(res)
  
  # Add gene symbol from our mapping
  res <- res %>%
    left_join(probe_map, by = "PROBEID") %>%
    select(PROBEID, SYMBOL, logFC, AveExpr,
           t, P.Value, adj.P.Val, B) %>%
    arrange(adj.P.Val)
  
  # Store in list
  all_results[[contrast_name]] <- res
  
  # Save to CSV
  filename <- paste0(
    "results/differential_expression/DE_",
    contrast_name, ".csv"
  )
  write.csv(res, filename, row.names = FALSE)
  
  cat("✓ Saved:", contrast_name,
      "(", nrow(res), "probes )\n")
}




# STEP 13: EXTRACT GLYCOGEN GENE RESULTS

# Our key genes grouped by biological function
glycogen_genes <- c(
  # Synthesis
  "GYS1", "GYS2", "UGP2", "GBE1", "PPP1R3C", "PPP1R3D",
  # Degradation
  "PYGB", "PYGL",
  # Regulation
  "GSK3B", "AKT1",
  # Warburg
  "HK2", "LDHA", "PKM", "G6PD", "SLC2A1",
  # Differentiation
  "ESR1", "KRT4", "KRT13"
)

# Extract from Cancer vs Normal (our primary comparison)
cancer_results <- all_results[["Cancer_vs_Normal"]]

glycogen_results <- cancer_results %>%
  filter(SYMBOL %in% glycogen_genes) %>%
  filter(!is.na(SYMBOL)) %>%
  # Keep best probe per gene (lowest p-value)
  group_by(SYMBOL) %>%
  slice_min(P.Value, n = 1) %>%
  ungroup() %>%
  arrange(SYMBOL) %>%
  mutate(
    logFC     = round(logFC,     3),
    adj.P.Val = round(adj.P.Val, 4),
    Direction = case_when(
      logFC > 0 ~ "UP in Cancer",
      logFC < 0 ~ "DOWN in Cancer",
      TRUE      ~ "No change"
    ),
    Significant = ifelse(adj.P.Val < 0.05, "YES", "NO")
  ) %>%
  dplyr::select(SYMBOL, logFC, adj.P.Val,
                Direction, Significant)

cat("Glycogen gene results — Cancer vs Normal:\n\n")
print(glycogen_results, n = Inf)

# Save
write.csv(
  glycogen_results,
  "results/differential_expression/glycogen_Cancer_vs_Normal.csv",
  row.names = FALSE
)


# VOLCANO PLOT — CANCER vs NORMAL

library(ggrepel)

cat("Building volcano plot...\n\n")

# Prepare full dataset for plotting
volcano_data <- all_results[["Cancer_vs_Normal"]] %>%
  filter(!is.na(SYMBOL)) %>%
  # Keep one row per gene symbol
  group_by(SYMBOL) %>%
  slice_min(P.Value, n = 1) %>%
  ungroup() %>%
  mutate(
    neg_log10_p = -log10(adj.P.Val),
    
    # Classify each gene for colouring
    gene_class = case_when(
      SYMBOL %in% c("GYS1","GYS2","GBE1",
                    "UGP2","PPP1R3C","PPP1R3D") ~ "Glycogen Synthesis",
      SYMBOL %in% c("PYGB","PYGL")              ~ "Glycogen Degradation",
      SYMBOL %in% c("HK2","LDHA","PKM",
                    "G6PD","SLC2A1")            ~ "Warburg Pathway",
      SYMBOL %in% c("ESR1","KRT4","KRT13")      ~ "Differentiation",
      adj.P.Val < 0.05 & abs(logFC) > 1         ~ "Other Significant",
      TRUE                                       ~ "Background"
    ),
    
    # Label only our key glycogen genes
    label = case_when(
      SYMBOL %in% glycogen_genes ~ SYMBOL,
      TRUE ~ NA_character_
    )
  )

# Colours
volcano_colours <- c(
  "Glycogen Synthesis"   = "#27AE60",
  "Glycogen Degradation" = "#1E8449",
  "Warburg Pathway"      = "#E74C3C",
  "Differentiation"      = "#2980B9",
  "Other Significant"    = "#AAB7B8",
  "Background"           = "#E8E8E8"
)

# Point sizes
volcano_sizes <- c(
  "Glycogen Synthesis"   = 4,
  "Glycogen Degradation" = 4,
  "Warburg Pathway"      = 4,
  "Differentiation"      = 4,
  "Other Significant"    = 1,
  "Background"           = 0.4
)

# Build plot
volcano_plot <- ggplot(
  volcano_data,
  aes(x = logFC, y = neg_log10_p,
      colour = gene_class,
      size   = gene_class)
) +
  
  # All background genes first
  geom_point(alpha = 0.5) +
  
  # Label our key genes
  geom_label_repel(
    aes(label = label),
    size         = 3,
    fontface     = "bold",
    fill         = "white",
    colour       = "black",
    alpha        = 0.85,
    box.padding  = 0.5,
    max.overlaps = 30,
    na.rm        = TRUE,
    segment.size = 0.3
  ) +
  
  # Reference lines
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             colour = "grey40", size = 0.4) +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed",
             colour = "grey40", size = 0.4) +
  
  # Colours and sizes
  scale_colour_manual(values = volcano_colours) +
  scale_size_manual(values = volcano_sizes) +
  
  # Labels
  labs(
    title   = "Differential Expression: Cancer vs Normal Cervical Epithelium",
    subtitle = paste0(
      "Green = glycogen synthesis  |  Red = Warburg  |  ",
      "Blue = differentiation\n",
      "Dashed lines: FDR = 0.05 and |logFC| = 1"
    ),
    x       = "Log2 Fold Change (Cancer vs Normal)",
    y       = "-Log10 Adjusted P-value",
    colour  = "Gene Category",
    caption = "GSE63514  |  MolVIA Program — Publication 1"
  ) +
  
  theme_classic() +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9, colour = "grey40"),
    plot.caption    = element_text(size = 8, colour = "grey60"),
    legend.position = "right",
    legend.title    = element_text(face = "bold")
  ) +
  guides(size = "none")

# Save
ggsave(
  "figures/volcano_plots/volcano_Cancer_vs_Normal.png",
  volcano_plot,
  width  = 12,
  height = 8,
  dpi    = 300
)

print(volcano_plot)
