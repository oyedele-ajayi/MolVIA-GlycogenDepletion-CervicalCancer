# ============================================================
# STEP 1b: CREATE PROJECT FOLDER STRUCTURE
# ============================================================
# We are creating folders to organise our work neatly
# Each folder has one clear purpose

# Create main folders
dir.create("data",        showWarnings = FALSE)
dir.create("results",     showWarnings = FALSE)
dir.create("figures",     showWarnings = FALSE)
dir.create("scripts",     showWarnings = FALSE)
dir.create("reports",     showWarnings = FALSE)

# Create sub-folders inside data
dir.create("data/raw",       showWarnings = FALSE)
dir.create("data/processed", showWarnings = FALSE)
dir.create("data/gene_sets", showWarnings = FALSE)

# Create sub-folders inside results
dir.create("results/differential_expression", showWarnings = FALSE)
dir.create("results/pathway_analysis",        showWarnings = FALSE)
dir.create("results/GDI_scores",              showWarnings = FALSE)

# Create sub-folders inside figures
dir.create("figures/quality_control", showWarnings = FALSE)
dir.create("figures/heatmaps",        showWarnings = FALSE)
dir.create("figures/volcano_plots",   showWarnings = FALSE)
dir.create("figures/GDI_plots",       showWarnings = FALSE)

# Confirm all folders were created
cat("✓ Folder structure created successfully!\n\n")
cat("Your project folders:\n")
list.dirs(recursive = TRUE)


# ============================================================
# STEP 2a: INSTALL BIOCMANAGER
# ============================================================
# BiocManager lets us install Bioconductor packages
# Bioconductor = specialised R packages for biology

# Check if BiocManager is already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
  cat("✓ BiocManager installed\n")
} else {
  cat("✓ BiocManager already installed\n")
}

# Check BiocManager version
BiocManager::version()


# ============================================================
# STEP 2b: INSTALL ALL PACKAGES
# ============================================================

# First — standard R packages from CRAN
cran_packages <- c(
  "tidyverse",    # Data manipulation and plotting
  "ggplot2",      # Making beautiful figures
  "ggrepel",      # Labels on figures that don't overlap
  "pheatmap",     # Heatmap figures
  "RColorBrewer", # Colour palettes for figures
  "openxlsx",     # Save results to Excel files
  "here",         # Helps find files in our project
  "janitor",      # Clean up data tables
  "corrplot"      # Correlation plots
)

# Install CRAN packages
cat("Installing standard packages...\n")
install.packages(cran_packages)
cat("✓ Standard packages installed\n\n")

# Second — Bioconductor packages (biology-specific)
bioc_packages <- c(
  "GEOquery",       # Download data from GEO database
  "Biobase",        # Core biology data structures
  "limma",          # Differential expression analysis
  "affy",           # Affymetrix microarray processing
  "oligo",          # Alternative microarray processing
  "hgu133plus2.db", # Annotation for our microarray type
  "annotate",       # Gene annotation tools
  "GSVA",           # Gene set variation analysis (for GDI)
  "GSEABase"        # Gene set handling
)

# Install Bioconductor packages
cat("Installing Bioconductor packages...\n")
BiocManager::install(bioc_packages, ask = FALSE)
cat("✓ Bioconductor packages installed\n\n")
