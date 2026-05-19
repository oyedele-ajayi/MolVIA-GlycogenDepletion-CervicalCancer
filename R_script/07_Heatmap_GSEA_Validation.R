# SESSION 4

rm(list = ls())

# Load all libraries needed
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(hgu133plus2.db)
library(AnnotationDbi)


# LOAD RAW DATA

load("data/raw/GSE63514_downloaded.RData")

# Extract core objects
gse_data    <- gse63514[[1]]
sample_info <- pData(gse_data)
expr_matrix <- exprs(gse_data)

# Rename columns
colnames(expr_matrix) <- sample_info$title

cat("✓ Raw data loaded\n")
cat("  Probes :", nrow(expr_matrix), "\n")
cat("  Samples:", ncol(expr_matrix), "\n\n")



# REBUILD SAMPLE ANNOTATION

cat("Building sample annotation...\n")

correct_order <- c("Normal","CIN1","CIN2","CIN3","Cancer")

sample_annotation <- data.frame(
  sample_id = rownames(sample_info),
  title     = sample_info$title,
  group     = factor(
    str_remove(sample_info$title, "-\\d+$"),
    levels  = correct_order,
    ordered = TRUE
  ),
  stringsAsFactors = FALSE
)

cat("✓ Sample annotation built\n")
print(table(sample_annotation$group))
cat("\n")



# REBUILD PROBE MAP

cat("Building probe map...\n")

probe_map <- AnnotationDbi::select(
  hgu133plus2.db,
  keys    = rownames(expr_matrix),
  columns = c("PROBEID","SYMBOL"),
  keytype = "PROBEID"
) %>%
  filter(!is.na(SYMBOL))

cat("✓ Probe map built\n")
cat("  Mapped probes:", nrow(probe_map), "\n")
cat("  Unique genes :", length(unique(probe_map$SYMBOL)), "\n\n")



# BUILD GENE-LEVEL EXPRESSION MATRIX

expr_gene_level <- expr_matrix %>%
  as.data.frame() %>%
  mutate(PROBEID = rownames(.)) %>%
  left_join(probe_map, by = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  dplyr::select(-PROBEID) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("SYMBOL")

cat("  Unique genes:", nrow(expr_gene_level), "\n")
cat("  Samples     :", ncol(expr_gene_level), "\n\n")



# BUILD HEATMAP OBJECTS

# Genes for heatmap
heatmap_genes <- c(
  "GYS1","GYS2","GBE1","UGP2","PPP1R3C","PPP1R3D",
  "PYGB","PYGL",
  "GSK3B","AKT1",
  "HK2","LDHA","PKM","SLC2A1","G6PD",
  "ESR1","KRT4","KRT13"
)

# Extract and order
heatmap_matrix <- expr_gene_level[
  rownames(expr_gene_level) %in% heatmap_genes, ]

col_order <- order(sample_annotation$group)
heatmap_matrix <- heatmap_matrix[, col_order]
sample_annotation_ordered <- sample_annotation[col_order, ]

# Scale rows
heatmap_scaled <- t(scale(t(heatmap_matrix)))

cat("✓ Heatmap matrix built\n")
cat("  Genes  :", nrow(heatmap_scaled), "\n")
cat("  Samples:", ncol(heatmap_scaled), "\n\n")

# Column annotation
col_annotation <- data.frame(
  Grade     = sample_annotation_ordered$group,
  row.names = sample_annotation_ordered$title
)

# Row annotation
row_annotation <- data.frame(
  Function = case_when(
    rownames(heatmap_scaled) %in%
      c("GYS1","GYS2","GBE1","UGP2",
        "PPP1R3C","PPP1R3D")       ~ "Glycogen Synthesis",
    rownames(heatmap_scaled) %in%
      c("PYGB","PYGL")             ~ "Glycogen Degradation",
    rownames(heatmap_scaled) %in%
      c("GSK3B","AKT1")            ~ "Regulation",
    rownames(heatmap_scaled) %in%
      c("HK2","LDHA","PKM",
        "SLC2A1","G6PD")           ~ "Warburg Pathway",
    rownames(heatmap_scaled) %in%
      c("ESR1","KRT4","KRT13")     ~ "Differentiation",
    TRUE ~ "Other"
  ),
  row.names = rownames(heatmap_scaled)
)

# Annotation colours
annotation_colours <- list(
  Grade = c(
    "Normal" = "#2ECC71",
    "CIN1"   = "#F39C12",
    "CIN2"   = "#E67E22",
    "CIN3"   = "#E74C3C",
    "Cancer" = "#8E44AD"
  ),
  Function = c(
    "Glycogen Synthesis"   = "#1A9850",
    "Glycogen Degradation" = "#91CF60",
    "Regulation"           = "#FEE08B",
    "Warburg Pathway"      = "#D73027",
    "Differentiation"      = "#4575B4"
  )
)

# Group colours for later steps
group_colours <- c(
  "Normal" = "#2ECC71",
  "CIN1"   = "#F39C12",
  "CIN2"   = "#E67E22",
  "CIN3"   = "#E74C3C",
  "Cancer" = "#8E44AD"
)



# DRAW IMPROVED HEATMAP

heatmap_colours_v2 <- colorRampPalette(
  rev(brewer.pal(11, "RdBu"))
)(100)

dir.create("figures/heatmaps", showWarnings = FALSE,
           recursive = TRUE)

png(
  filename = "figures/heatmaps/glycogen_heatmap_publication.png",
  width    = 16,
  height   = 10,
  units    = "in",
  res      = 300
)

pheatmap(
  mat               = heatmap_scaled,
  color             = heatmap_colours_v2,
  cluster_cols      = FALSE,
  show_colnames     = FALSE,
  cluster_rows      = TRUE,
  show_rownames     = TRUE,
  fontsize_row      = 12,
  fontsize          = 11,
  annotation_col    = col_annotation,
  annotation_row    = row_annotation,
  annotation_colors = annotation_colours,
  border_color      = NA,
  breaks            = seq(-2.5, 2.5, length.out = 101),
  main              = paste0(
    "Glycogen Metabolic Reprogramming Across Cervical Disease Grades\n",
    "GSE63514 | n=128 | Row z-score | ",
    "Normal → CIN1 → CIN2 → CIN3 → Cancer"
  ),
  gaps_col = cumsum(c(24, 14, 22, 40)),
  legend   = TRUE
)

dev.off()

cat("✓ Heatmap saved:\n")
cat("  figures/heatmaps/glycogen_heatmap_publication.png\n\n")




# SAVE COMPLETE SESSION 4 WORKSPACE

save(
  expr_matrix,
  expr_gene_level,
  sample_annotation,
  sample_annotation_ordered,
  probe_map,
  correct_order,
  group_colours,
  heatmap_scaled,
  heatmap_matrix,
  col_annotation,
  row_annotation,
  annotation_colours,
  file = "data/processed/session04_workspace.RData"
)

cat("  data/processed/session04_workspace.RData\n\n")

cat("=== REBUILD COMPLETE ===\n\n")
cat("Objects now available:\n")
cat("  ✓ expr_matrix\n")
cat("  ✓ expr_gene_level\n")
cat("  ✓ sample_annotation\n")
cat("  ✓ probe_map\n")
cat("  ✓ heatmap_scaled\n")
cat("  ✓ group_colours\n")



# LOAD fgsea PACKAGE


# Install if needed
if (!requireNamespace("fgsea", quietly = TRUE)) {
  BiocManager::install("fgsea", ask = FALSE)
}

library(fgsea)

cat("fgsea version:", as.character(packageVersion("fgsea")), "\n")



# BUILD RANKED GENE LIST

# Load session 3 results
# We need all_results from limma analysis
# Load from the saved CSV since all_results
# may not be in session04_workspace

cancer_de <- read.csv(
  "results/differential_expression/DE_Cancer_vs_Normal.csv",
  stringsAsFactors = FALSE
)

# Remove rows with no gene symbol
cancer_de <- cancer_de %>%
  filter(!is.na(SYMBOL)) %>%
  filter(SYMBOL != "")

cat("After removing missing symbols:", nrow(cancer_de), "rows\n\n")

# Calculate ranking metric
# signed -log10(p) preserves direction and magnitude
cancer_de <- cancer_de %>%
  mutate(
    rank_metric = sign(logFC) * (-log10(P.Value))
  )

# For genes with multiple probes — keep best ranked
cancer_ranked <- cancer_de %>%
  group_by(SYMBOL) %>%
  slice_max(abs(rank_metric), n = 1) %>%
  ungroup() %>%
  arrange(desc(rank_metric))

# Build named numeric vector (required by fgsea)
ranked_list <- setNames(
  cancer_ranked$rank_metric,
  cancer_ranked$SYMBOL
)

cat("=== RANKED GENE LIST SUMMARY ===\n\n")
cat("Total genes ranked     :", length(ranked_list), "\n")
cat("Most upregulated gene  :", names(ranked_list)[1],
    "(score:", round(ranked_list[1], 2), ")\n")
cat("Most downregulated gene:", names(ranked_list)[length(ranked_list)],
    "(score:", round(ranked_list[length(ranked_list)], 2), ")\n\n")

# Where do our key genes rank?
cat("Position of key glycogen genes in ranking:\n\n")

key_check <- c("PPP1R3C","GYS2","GYS1",
               "SLC2A1","HK2","ESR1",
               "KRT4","KRT13")

for (gene in key_check) {
  if (gene %in% names(ranked_list)) {
    position <- which(names(ranked_list) == gene)
    score    <- round(ranked_list[gene], 3)
    total    <- length(ranked_list)
    pct      <- round(position/total * 100, 1)
    cat(sprintf("  %-10s rank %5d of %d (%5.1f%%)  score: %7.3f\n",
                gene, position, total, pct, score))
  } else {
    cat(sprintf("  %-10s not found in ranked list\n", gene))
  }
}



# STEP RUN GSEA


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

cat("✓ glycogen_modules defined\n")
for (mod in names(glycogen_modules)) {
  cat(sprintf("  %-25s: %d genes\n",
              mod, length(glycogen_modules[[mod]])))
}
cat("\n")

# Use fgseaMultilevel (no nperm argument)

set.seed(42)

cat("Running fgseaMultilevel...\n")

gsea_results <- fgseaMultilevel(
  pathways = glycogen_modules,
  stats    = ranked_list,
  minSize  = 5,
  maxSize  = 500
)


# Clean results table
gsea_clean <- gsea_results %>%
  arrange(NES) %>%
  mutate(
    NES       = round(NES,  3),
    pval      = round(pval, 5),
    padj      = round(padj, 5),
    Direction = ifelse(
      NES < 0,
      "Downregulated",
      "Upregulated"
    )
  ) %>%
  dplyr::select(pathway, NES, pval, padj, Direction)

print(gsea_clean)

# Significance summary
for (i in 1:nrow(gsea_clean)) {
  sig <- ifelse(gsea_clean$padj[i] < 0.05, "✓ SIG", "— ns")
  cat(sprintf("  %-25s NES = %6.3f  padj = %.5f  %s\n",
              gsea_clean$pathway[i],
              gsea_clean$NES[i],
              gsea_clean$padj[i],
              sig))
}

# Save
write.csv(
  gsea_clean,
  "results/pathway_analysis/GSEA_glycogen_Cancer_vs_Normal.csv",
  row.names = FALSE
)

cat("\n✓ GSEA results saved to:\n")
cat("  results/pathway_analysis/GSEA_glycogen_Cancer_vs_Normal.csv\n")




# GSEA ENRICHMENT PLOTS

library(ggplot2)

cat("Creating GSEA enrichment plots...\n\n")

dir.create("figures/pathway_analysis",
           showWarnings = FALSE,
           recursive    = TRUE)

# ── Plot 1: Differentiation (significant) ──────────────────

p1 <- plotEnrichment(
  pathway = glycogen_modules[["Differentiation"]],
  stats   = ranked_list
) +
  labs(
    title    = "GSEA — Differentiation Gene Set",
    subtitle = paste0(
      "NES = -1.876  |  padj = 0.003  |  ",
      "Significantly depleted in Cancer vs Normal"
    ),
    x = "Rank in gene list (Cancer vs Normal)",
    y = "Enrichment Score",
    caption = "GSE63514 | MolVIA Program — Publication 1"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    plot.caption  = element_text(size = 8, colour = "grey60")
  )

ggsave(
  "figures/pathway_analysis/GSEA_Differentiation.png",
  p1, width = 8, height = 5, dpi = 300
)

print(p1)
cat("✓ Differentiation GSEA plot saved\n\n")

# ── Plot 2: Glycogen Synthesis (trending) ──────────────────

p2 <- plotEnrichment(
  pathway = glycogen_modules[["Glycogen_Synthesis"]],
  stats   = ranked_list
) +
  labs(
    title    = "GSEA — Glycogen Synthesis Gene Set",
    subtitle = paste0(
      "NES = -1.477  |  padj = 0.114  |  ",
      "Negative enrichment trend in Cancer vs Normal"
    ),
    x = "Rank in gene list (Cancer vs Normal)",
    y = "Enrichment Score",
    caption = "GSE63514 | MolVIA Program — Publication 1"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    plot.caption  = element_text(size = 8, colour = "grey60")
  )

ggsave(
  "figures/pathway_analysis/GSEA_Glycogen_Synthesis.png",
  p2, width = 8, height = 5, dpi = 300
)

print(p2)


# ── Plot 3: Warburg Pathway ────

p3 <- plotEnrichment(
  pathway = glycogen_modules[["Warburg_Pathway"]],
  stats   = ranked_list
) +
  labs(
    title    = "GSEA — Warburg Pathway Gene Set",
    subtitle = paste0(
      "NES = +1.261  |  padj = 0.278  |  ",
      "Positive enrichment trend in Cancer vs Normal"
    ),
    x = "Rank in gene list (Cancer vs Normal)",
    y = "Enrichment Score",
    caption = "GSE63514 | MolVIA Program — Publication 1"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    plot.caption  = element_text(size = 8, colour = "grey60")
  )

ggsave(
  "figures/pathway_analysis/GSEA_Warburg_Pathway.png",
  p3, width = 8, height = 5, dpi = 300
)

print(p3)
cat("✓ Warburg GSEA plot saved\n\n")
cat("All GSEA plots saved to figures/pathway_analysis/\n")





# DOWNLOAD VALIDATION DATASET GSE9750


library(GEOquery)

gse9750 <- getGEO(
  GEO       = "GSE9750",
  destdir   = "data/raw",
  GSEMatrix = TRUE,
  AnnotGPL  = TRUE
)

# Extract data
gse9750_data    <- gse9750[[1]]
sample_info_v   <- pData(gse9750_data)
expr_matrix_v   <- exprs(gse9750_data)

cat("=== GSE9750 OVERVIEW ===\n\n")
cat("Samples:", ncol(expr_matrix_v), "\n")
cat("Probes :", nrow(expr_matrix_v), "\n\n")

# sample titles
cat("Sample titles (first 10):\n")
print(sample_info_v$title[1:10])

cat("\nSource names:\n")
print(table(sample_info_v$source_name_ch1))




# CLEAN AND ANNOTATE GSE9750

# Rename columns to sample titles
colnames(expr_matrix_v) <- sample_info_v$title

# Classify each sample
sample_annotation_v <- data.frame(
  title  = sample_info_v$title,
  source = sample_info_v$source_name_ch1,
  group  = case_when(
    grepl("Normal", sample_info_v$source_name_ch1,
          ignore.case = TRUE)       ~ "Normal",
    grepl("cell line", sample_info_v$source_name_ch1,
          ignore.case = TRUE)       ~ "Cell_line",
    grepl("cancer", sample_info_v$source_name_ch1,
          ignore.case = TRUE)       ~ "Cancer_tissue",
    TRUE ~ "Unknown"
  ),
  stringsAsFactors = FALSE
)

print(table(sample_annotation_v$group))

# Exclude cell lines — keep only primary tissue
sample_annotation_v <- sample_annotation_v %>%
  filter(group != "Cell_line")

cat("\nAfter excluding cell lines:\n\n")
print(table(sample_annotation_v$group))

# Filter expression matrix to keep only primary tissue
expr_matrix_v <- expr_matrix_v[,
                               sample_annotation_v$title]

cat("\nValidation matrix dimensions:\n")
cat("  Probes :", nrow(expr_matrix_v), "\n")
cat("  Samples:", ncol(expr_matrix_v), "\n")




# FIX UNKNOWN SAMPLE AND FINALISE

# Show the unknown sample
unknown_sample <- sample_annotation_v %>%
  filter(group == "Unknown")

print(unknown_sample)

# reclassify as Cancer_tissue
# "Cervical cacner_654T" is clearly a cancer sample
# with a typo in GEO annotation
sample_annotation_v <- sample_annotation_v %>%
  mutate(
    group = ifelse(group == "Unknown",
                   "Cancer_tissue",
                   group)
  )

print(table(sample_annotation_v$group))

# Set as factor with correct order
sample_annotation_v$group <- factor(
  sample_annotation_v$group,
  levels = c("Normal", "Cancer_tissue")
)

cat("\n✓ Groups finalised\n")
cat("  Normal        :", sum(sample_annotation_v$group == "Normal"), "\n")
cat("  Cancer_tissue :", sum(sample_annotation_v$group == "Cancer_tissue"), "\n")
cat("  Total         :", nrow(sample_annotation_v), "\n")




# BUILD GENE-LEVEL MATRIX FOR GSE9750


# Install GPL96 annotation if needed
if (!requireNamespace("hgu133a.db", quietly = TRUE)) {
  BiocManager::install("hgu133a.db", ask = FALSE)
}

library(hgu133a.db)

cat("Mapping probes to gene symbols (GPL96 platform)...\n\n")

# Map probes to gene symbols
probe_map_v <- AnnotationDbi::select(
  hgu133a.db,
  keys    = rownames(expr_matrix_v),
  columns = c("PROBEID", "SYMBOL"),
  keytype = "PROBEID"
) %>%
  filter(!is.na(SYMBOL))

cat("  Mapped probes:", nrow(probe_map_v), "\n")
cat("  Unique genes :", length(unique(probe_map_v$SYMBOL)), "\n\n")

# Collapse probes to gene level
expr_gene_level_v <- expr_matrix_v %>%
  as.data.frame() %>%
  mutate(PROBEID = rownames(.)) %>%
  left_join(probe_map_v, by = "PROBEID") %>%
  filter(!is.na(SYMBOL)) %>%
  dplyr::select(-PROBEID) %>%
  group_by(SYMBOL) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("SYMBOL")

cat("  Unique genes:", nrow(expr_gene_level_v), "\n")
cat("  Samples     :", ncol(expr_gene_level_v), "\n\n")

# Check if key glycogen genes are present
key_genes_check <- c("GYS1","GYS2","PPP1R3C",
                     "SLC2A1","HK2","ESR1",
                     "KRT4","KRT13")

for (gene in key_genes_check) {
  present <- gene %in% rownames(expr_gene_level_v)
  symbol  <- ifelse(present, "✓", "✗")
  cat(" ", symbol, gene, "\n")
}





# ssGSEA + GDI IN VALIDATION DATASET

# LOAD GSVA 
library(GSVA)

# Define the same modules (must be identical)
glycogen_modules <- list(
  Glycogen_Synthesis = c(
    "GYS1","GYS2","GBE1","UGP2",
    "PGM1","PPP1R3A","PPP1R3B",
    "PPP1R3C","PPP1R3D"
  ),
  Glycogen_Degradation = c(
    "PYGB","PYGL","PYGM","AGL",
    "PHKA1","PHKB","PHKG1"
  ),
  Warburg_Pathway = c(
    "HK2","LDHA","PKM","SLC2A1",
    "SLC2A3","PFKL","ENO1","G6PD"
  ),
  Differentiation = c(
    "ESR1","ESR2","TP63","KRT4",
    "KRT13","KRT14","IVL","SPRR1A","SPRR1B"
  )
)

# Run ssGSEA (version-safe)
set.seed(42)

if (packageVersion("GSVA") >= "1.44") {
  
  ssgsea_param_v <- ssgseaParam(
    exprData = as.matrix(expr_gene_level_v),
    geneSets = glycogen_modules
  )
  ssgsea_scores_v <- gsva(ssgsea_param_v, verbose = FALSE)
  
} else {
  
  ssgsea_scores_v <- gsva(
    expr          = as.matrix(expr_gene_level_v),
    gset.idx.list = glycogen_modules,
    method        = "ssgsea",
    kcdf          = "Gaussian",
    verbose       = FALSE
  )
}

cat("Score matrix:", nrow(ssgsea_scores_v), "modules x",
    ncol(ssgsea_scores_v), "samples\n\n")

# Compute GDI
synthesis_v <- ssgsea_scores_v["Glycogen_Synthesis", ]
degrad_v    <- ssgsea_scores_v["Glycogen_Degradation", ]
warburg_v   <- ssgsea_scores_v["Warburg_Pathway", ]
diff_v      <- ssgsea_scores_v["Differentiation", ]

GDI_v <- (warburg_v + degrad_v) / (synthesis_v + diff_v)

# Build table
gdi_table_v <- data.frame(
  sample = names(GDI_v),
  GDI    = as.numeric(GDI_v),
  group  = sample_annotation_v$group
)

# Summary

gdi_summary_v <- gdi_table_v %>%
  group_by(group) %>%
  summarise(
    n          = n(),
    mean_GDI   = round(mean(GDI), 4),
    median_GDI = round(median(GDI), 4),
    sd_GDI     = round(sd(GDI), 4),
    .groups = "drop"
  )

print(gdi_summary_v)

# Quick test: Cancer vs Normal
cat("\nWilcoxon test (Cancer vs Normal):\n")
wt_v <- wilcox.test(GDI ~ group, data = gdi_table_v)
cat("  p-value:", format(wt_v$p.value, scientific = TRUE), "\n")

if (wt_v$p.value < 0.05) {
  cat("  ✓ SIGNIFICANT — GDI higher in cancer in validation cohort\n")
} else {
  cat(" Not significant — report to PI\n")
}




# VALIDATION ROC CURVE AND AUC


library(pROC)

# Binary outcome: Cancer=1, Normal=0
gdi_table_v <- gdi_table_v %>%
  mutate(
    cancer_binary = ifelse(group == "Cancer_tissue", 1, 0)
  )

# ROC curve
roc_validation <- roc(
  gdi_table_v$cancer_binary,
  gdi_table_v$GDI,
  quiet = TRUE
)

auc_validation <- round(auc(roc_validation), 3)
ci_validation  <- round(ci.auc(roc_validation), 3)

cat("Validation dataset (GSE9750):\n")
cat(sprintf("  AUC = %s (95%% CI: %s - %s)\n",
            auc_validation,
            ci_validation[1],
            ci_validation[3]))

# Interpret
cat("\n  Interpretation: GDI for Cancer detection → AUC =",
    auc_validation)
if (auc_validation >= 0.80) {
  cat("  ✓ EXCELLENT discrimination\n")
} else if (auc_validation >= 0.70) {
  cat("  ✓ GOOD discrimination\n")
} else if (auc_validation >= 0.60) {
  cat("  ~ MODERATE discrimination\n")
} else {
  cat("POOR")
}

# Cross-dataset AUC comparison
cat("\n=== CROSS-DATASET AUC COMPARISON ===\n\n")
cat(sprintf("  %-20s AUC = %s (95%% CI: %s - %s)\n",
            "GSE63514 CIN3+:",
            0.814, 0.741, 0.888))
cat(sprintf("  %-20s AUC = %s (95%% CI: %s - %s)\n",
            "GSE9750 Cancer:",
            auc_validation,
            ci_validation[1],
            ci_validation[3]))




# VALIDATION GDI BOXPLOT

cat("\nCreating validation GDI plot...\n\n")

# Validation colours
validation_colours <- c(
  "Normal"        = "#2ECC71",
  "Cancer_tissue" = "#8E44AD"
)

# Validation boxplot
gdi_validation_plot <- ggplot(
  gdi_table_v,
  aes(x = group, y = GDI, fill = group)
) +
  geom_boxplot(
    alpha = 0.7,
    width = 0.5,
    outlier.shape = 21,
    outlier.size  = 2
  ) +
  geom_jitter(
    width  = 0.12,
    size   = 1.8,
    alpha  = 0.5,
    colour = "black"
  ) +
  stat_summary(
    fun    = median,
    geom   = "point",
    shape  = 23,
    size   = 4.5,
    fill   = "white",
    colour = "black"
  ) +
  scale_fill_manual(values = validation_colours) +
  scale_x_discrete(
    labels = c(
      "Normal"        = "Normal\n(n=24)",
      "Cancer_tissue" = "Cancer\n(n=33)"
    )
  ) +
  annotate(
    "text",
    x = 1.5, y = max(gdi_table_v$GDI) * 0.85,
    label = paste0(
      "AUC = ", auc_validation,
      "\np = ", format(wt_v$p.value,
                       scientific = TRUE,
                       digits = 2)
    ),
    size     = 4,
    fontface = "bold",
    colour   = "grey20"
  ) +
  labs(
    title    = "Validation: GDI in Independent Cohort (GSE9750)",
    subtitle = paste0(
      "Normal (n=24) vs Primary Cervical Cancer (n=33) | ",
      "Independent patient cohort | GPL96 platform"
    ),
    x       = "Sample Group",
    y       = "Glycogen Depletion Index",
    caption = "GSE9750 | MolVIA Program — Publication 1 | Validation"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey40"),
    plot.caption  = element_text(size = 8, colour = "grey60"),
    axis.title    = element_text(face = "bold"),
    legend.position = "none"
  )

ggsave(
  "figures/GDI_plots/GDI_validation_GSE9750.png",
  gdi_validation_plot,
  width  = 7,
  height = 7,
  dpi    = 300
)

print(gdi_validation_plot)




# STEP 21h: COMBINED DISCOVERY + VALIDATION FIGURE

library(patchwork)

# Reload discovery GDI plot data
# (from session03 results)
load("data/processed/session03_complete.RData")

# Rebuild discovery plot (zoomed)
y_upper_disc <- quantile(gdi_table$GDI, 0.95,
                         na.rm = TRUE)

discovery_plot <- ggplot(
  gdi_table,
  aes(x = group, y = GDI, fill = group)
) +
  geom_boxplot(
    alpha = 0.7, width = 0.5,
    outlier.shape = 21, outlier.size = 1.5
  ) +
  geom_jitter(
    width = 0.1, size = 1.5,
    alpha = 0.4, colour = "black"
  ) +
  stat_summary(
    fun = median, geom = "point",
    shape = 23, size = 4,
    fill = "white", colour = "black"
  ) +
  scale_fill_manual(values = group_colours) +
  coord_cartesian(ylim = c(0, y_upper_disc)) +
  labs(
    title    = "A. Discovery Cohort (GSE63514)",
    subtitle = paste0(
      "n=128 | Kruskal-Wallis p=2.19×10⁻¹⁰\n",
      "Spearman rho=0.621 | AUC CIN3+=0.814"
    ),
    x = "Disease Grade",
    y = "Glycogen Depletion Index"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    axis.title    = element_text(face = "bold"),
    legend.position = "none"
  )

# Validation plot
validation_plot_clean <- ggplot(
  gdi_table_v,
  aes(x = group, y = GDI, fill = group)
) +
  geom_boxplot(
    alpha = 0.7, width = 0.5,
    outlier.shape = 21, outlier.size = 1.5
  ) +
  geom_jitter(
    width = 0.1, size = 1.5,
    alpha = 0.4, colour = "black"
  ) +
  stat_summary(
    fun = median, geom = "point",
    shape = 23, size = 4,
    fill = "white", colour = "black"
  ) +
  scale_fill_manual(values = validation_colours) +
  scale_x_discrete(
    labels = c(
      "Normal"        = "Normal",
      "Cancer_tissue" = "Cancer"
    )
  ) +
  annotate(
    "text", x = 1.5,
    y = max(gdi_table_v$GDI) * 0.80,
    label = paste0("AUC = ", auc_validation,
                   "\np = 3.93×10⁻⁸"),
    size = 3.5, fontface = "bold", colour = "grey20"
  ) +
  labs(
    title    = "B. Validation Cohort (GSE9750)",
    subtitle = paste0(
      "n=57 | Wilcoxon p=3.93×10⁻⁸\n",
      "AUC = ", auc_validation,
      " | Independent cohort + platform"
    ),
    x = "Sample Group",
    y = "Glycogen Depletion Index"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 8, colour = "grey40"),
    axis.title    = element_text(face = "bold"),
    legend.position = "none"
  )

# Combine panels
combined_validation <- discovery_plot + validation_plot_clean +
  plot_annotation(
    title   = "Glycogen Depletion Index — Discovery and Independent Validation",
    subtitle = paste0(
      "GDI consistently elevated in cervical dysplasia/cancer ",
      "across two independent cohorts and microarray platforms"
    ),
    caption = "MolVIA Program — Publication 1",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, colour = "grey40"),
      plot.caption  = element_text(size = 8, colour = "grey60")
    )
  )

ggsave(
  "figures/GDI_plots/GDI_discovery_validation_combined.png",
  combined_validation,
  width  = 14,
  height = 7,
  dpi    = 300
)

print(combined_validation)

library(conflicted)
library(tibble)

# Resolve conflicts
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)
conflicts_prefer(base::setdiff)
conflicts_prefer(base::intersect)
conflicts_prefer(base::union)

# Safe print function that works for both data frames and tibbles
safe_print <- function(x, ...) {
  print(tibble::as_tibble(x), n = Inf, ...)
}


# MASTER RESULTS SUMMARY TABLE


cat("=== BUILDING MASTER RESULTS SUMMARY TABLE ===\n\n")

# ── Section 1: Key Glycogen Gene Results ──────

gene_results_table <- data.frame(
  Section   = "Differential Expression",
  Finding   = c(
    "PPP1R3C (GYS activator)",
    "GYS2 (glycogen synthase)",
    "KRT13 (differentiation)",
    "KRT4  (differentiation)",
    "ESR1  (estrogen receptor)",
    "SLC2A1/GLUT1 (Warburg)",
    "PYGB  (glycogen degradation)",
    "GSK3B (regulation)"
  ),
  Statistic = c(
    "logFC = -5.12",
    "logFC = -3.74",
    "logFC = -4.67",
    "logFC = -3.44",
    "logFC = -0.64",
    "logFC = +1.76",
    "logFC = -0.29",
    "logFC = -0.26"
  ),
  FDR = c(
    "< 0.0001",
    "< 0.0001",
    "< 0.0001",
    "< 0.0001",
    "< 0.0001",
    "0.0008",
    "0.0006",
    "< 0.0001"
  ),
  Direction = c(
    "DOWN in Cancer",
    "DOWN in Cancer",
    "DOWN in Cancer",
    "DOWN in Cancer",
    "DOWN in Cancer",
    "UP in Cancer",
    "DOWN in Cancer",
    "DOWN in Cancer"
  ),
  stringsAsFactors = FALSE
)

# ── Section 2: GDI Performance ─────────────────────────────

gdi_results_table <- data.frame(
  Section   = "GDI v1.0 Performance",
  Finding   = c(
    "GDI — Kruskal-Wallis (all grades)",
    "GDI — Spearman correlation",
    "GDI — CIN3 vs Normal",
    "GDI — Cancer vs Normal",
    "GDI — AUC for CIN2+",
    "GDI — AUC for CIN3+"
  ),
  Statistic = c(
    "χ² = 51.04, df = 4",
    "rho = 0.621",
    "Wilcoxon",
    "Wilcoxon",
    "AUC = 0.803",
    "AUC = 0.814"
  ),
  FDR = c(
    "p = 2.19×10⁻¹⁰",
    "p = 5.17×10⁻¹⁵",
    "p = 6.30×10⁻⁵",
    "p = 6.87×10⁻⁹",
    "95% CI: 0.724-0.882",
    "95% CI: 0.741-0.888"
  ),
  Direction = c(
    "GDI increases with grade",
    "Positive monotonic trend",
    "GDI higher in CIN3",
    "GDI higher in Cancer",
    "Excellent discrimination",
    "Excellent discrimination"
  ),
  stringsAsFactors = FALSE
)

# ──Section 3: GSEA Results ────────────────────────────────

gsea_results_table <- data.frame(
  Section   = "GSEA Pathway Enrichment",
  Finding   = c(
    "Differentiation pathway",
    "Glycogen Synthesis pathway",
    "Warburg pathway",
    "Glycogen Degradation pathway"
  ),
  Statistic = c(
    "NES = -1.876",
    "NES = -1.477",
    "NES = +1.261",
    "NES = -0.567"
  ),
  FDR = c(
    "padj = 0.003",
    "padj = 0.114",
    "padj = 0.278",
    "padj = 0.945"
  ),
  Direction = c(
    "Significantly depleted",
    "Trend toward depletion",
    "Trend toward activation",
    "Not enriched"
  ),
  stringsAsFactors = FALSE
)

# ── Section 4: Validation Results ──────────────────────────

validation_results_table <- data.frame(
  Section   = "Independent Validation (GSE9750)",
  Finding   = c(
    "GDI Normal vs Cancer",
    "GDI AUC (Cancer detection)"
  ),
  Statistic = c(
    "Normal median=0.539, Cancer median=0.906",
    "AUC = 0.895"
  ),
  FDR = c(
    "p = 3.93×10⁻⁸",
    "95% CI: 0.795-0.995"
  ),
  Direction = c(
    "GDI higher in Cancer",
    "Excellent discrimination"
  ),
  stringsAsFactors = FALSE
)

# ── Combine all sections ──────

master_table <- bind_rows(
  gene_results_table,
  gdi_results_table,
  gsea_results_table,
  validation_results_table
)

cat("=== PUBLICATION 1 MASTER RESULTS TABLE ===\n\n")
safe_print(master_table)



# PRINT AND SAVE MASTER TABLE

# convert to data frame before printing
print(as.data.frame(master_table))

cat("\nTotal rows:", nrow(master_table), "\n\n")

# Save as CSV
write.csv(
  master_table,
  "results/Publication1_master_results_table.csv",
  row.names = FALSE
)

cat("✓ CSV saved\n")

# Save as Excel
library(openxlsx)

wb <- createWorkbook()
addWorksheet(wb, "Publication1_Results")

writeData(wb, "Publication1_Results", master_table)

# Header style
headerStyle <- createStyle(
  fontColour     = "#FFFFFF",
  fgFill         = "#2C3E50",
  halign         = "CENTER",
  fontName       = "Arial",
  fontSize       = 11,
  textDecoration = "bold",
  border         = "Bottom"
)

addStyle(
  wb, "Publication1_Results",
  headerStyle,
  rows      = 1,
  cols      = 1:ncol(master_table),
  gridExpand = TRUE
)

# Alternate row shading
for (row in seq(2, nrow(master_table) + 1, by = 2)) {
  addStyle(
    wb, "Publication1_Results",
    createStyle(fgFill = "#EBF5FB"),
    rows      = row,
    cols      = 1:ncol(master_table),
    gridExpand = TRUE
  )
}

# Auto column widths
setColWidths(
  wb, "Publication1_Results",
  cols   = 1:ncol(master_table),
  widths = "auto"
)

saveWorkbook(
  wb,
  "results/Publication1_master_results_table.xlsx",
  overwrite = TRUE
)

cat("✓ Excel saved\n\n")
cat("Files saved to results/ folder:\n")
cat("  ✓ Publication1_master_results_table.csv\n")
cat("  ✓ Publication1_master_results_table.xlsx\n")






# FINAL COMPLETE WORKSPACE SAVE


cat("=== SAVING FINAL PUBLICATION 1 WORKSPACE ===\n\n")

save(
  # ── Core data objects ──────────────────────────────────
  expr_matrix,
  expr_gene_level,
  expr_matrix_v,
  expr_gene_level_v,
  
  # ── Sample annotations ─────────────────────────────────
  sample_annotation,
  sample_annotation_ordered,
  sample_annotation_v,
  sample_info,
  sample_info_v,
  correct_order,
  
  # ── Probe mappings ─────────────────────────────────────
  probe_map,
  probe_map_v,
  
  # ── Gene sets and colours ──────────────────────────────
  glycogen_modules,
  group_colours,
  heatmap_genes,
  
  # ── Heatmap objects ────────────────────────────────────
  heatmap_scaled,
  heatmap_matrix,
  col_annotation,
  row_annotation,
  annotation_colours,
  
  # ── Limma analysis objects ─────────────────────────────
  design_matrix,
  contrast_matrix,
  fit1,
  fit2,
  
  # ── Differential expression results ───────────────────
  all_results,
  glycogen_results,
  
  # ── ssGSEA scores ──────────────────────────────────────
  ssgsea_scores,
  ssgsea_scores_v,
  
  # ── GDI results ────────────────────────────────────────
  gdi_table,
  gdi_summary,
  gdi_table_v,
  gdi_summary_v,
  
  # ── Statistical tests ──────────────────────────────────
  kw_test,
  spearman_test,
  wt_v,
  
  # ── ROC objects ────────────────────────────────────────
  roc_cin2,
  roc_cin3,
  roc_validation,
  auc_cin2,
  auc_cin3,
  auc_validation,
  ci_validation,
  
  # ── GSEA objects ───────────────────────────────────────
  gsea_results,
  gsea_clean,
  ranked_list,
  
  # ── Summary table ──────────────────────────────────────
  master_table,
  
  file = "data/processed/Publication1_FINAL_workspace.RData"
)

cat("✓ Final workspace saved\n\n")
