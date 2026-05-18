# RELOAD OUR SAVED DATA
load("data/raw/GSE63514_downloaded.RData")

# Quick confirmation dataset attribute
cat("Samples:", ncol(gse_data), "\n")
cat("Probes:", nrow(gse_data), "\n")


# VIEW ALL 128 SAMPLE TITLE
# Print all titles
all_titles <- sample_info$title

for (i in 1:length(all_titles)) {
  cat(i, ":", all_titles[i], "\n")
}


# STEP 6c: EXTRACT DISEASE GROUP FROM SAMPLE TITLE

library(tidyverse)

# Extract group name (everything before the "-" and number)
sample_groups <- sample_info %>%
  # Keep the title column
  select(title) %>%
  # Create a new column called 'group'
  # We take everything before the last hyphen+number
  mutate(
    group = str_remove(title, "-\\d+$")
  )

print(sample_groups)


# COUNT SAMPLES IN EACH DISEASE GROUP

group_counts <- sample_groups %>%
  count(group) %>%
  rename(
    "Disease Group" = group,
    "Number of Samples" = n
  )

print(group_counts)

cat("\nTotal samples:", sum(group_counts$`Number of Samples`), "\n")

# Also show as a simple table
print(table(sample_groups$group))



# CREATE CLEAN SAMPLE ANNOTATION TABLE

# Define the correct disease progression order
disease_order <- c("Normal", "CIN1", "CIN2", "CIN3", "Cancer")

# Build the sample annotation table
sample_annotation <- data.frame(
  sample_id = rownames(sample_info),
  title = sample_info$title,
  group = str_remove(sample_info$title, "-\\d+$"),
  stringsAsFactors = FALSE
)

# Set group as an ORDERED factor (Normal is lowest, Cancer highest)
sample_annotation$group <- factor(
  sample_annotation$group,
  levels = disease_order,
  ordered = TRUE
)

print(sample_annotation)

print(summary(sample_annotation$group))




# EXTRACT AND LABEL EXPRESSION MATRIX

# Extract the numeric expression values
expr_matrix <- exprs(gse_data)

cat("Number of probes (rows)   :", nrow(expr_matrix), "\n")
cat("Number of samples (columns):", ncol(expr_matrix), "\n")

# Rename columns to match sample titles
colnames(expr_matrix) <- sample_info$title

print(colnames(expr_matrix)[1:6])



# CHECK EXPRESSION VALUE RANGE

min_val  <- round(min(expr_matrix), 3)
max_val  <- round(max(expr_matrix), 3)
mean_val <- round(mean(expr_matrix), 3)

cat("Minimum:", min_val, "\n")
cat("Maximum:", max_val, "\n")
cat("Mean   :", mean_val, "\n\n")



# QC BOXPLOT

library(ggplot2)
library(tidyverse)

# Select every 3rd sample to keep the plot readable
sample_subset <- seq(1, ncol(expr_matrix), by = 3)

# Prepare data for ggplot
boxplot_data <- expr_matrix[, sample_subset] %>%
  as.data.frame() %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    group = factor(
      str_remove(sample, "-\\d+$"),
      levels = c("Normal","CIN1","CIN2","CIN3","Cancer")
    )
  )
group_colours <- c(
  "Normal" = "#2ECC71",
  "CIN1"   = "#F39C12",
  "CIN2"   = "#E67E22",
  "CIN3"   = "#E74C3C",
  "Cancer" = "#8E44AD"
)

# Create boxplot
qc_boxplot <- ggplot(boxplot_data, aes(x = sample, y = expression, fill = group)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
  scale_fill_manual(values = group_colours) +
  labs(
    title = "QC Boxplot — GSE63514 Expression Distribution",
    subtitle = "All boxes at similar height = good quality data",
    x = "Sample", y = "Expression (log2)", fill = "Group"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, size = 6, hjust = 1),
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "top"
  )

# Save and display
ggsave("figures/quality_control/QC_01_boxplot.png", qc_boxplot, width = 14, height = 6, dpi = 300)
print(qc_boxplot)




# QC DENSITY PLOT

density_data <- expr_matrix[, sample_subset] %>%
  as.data.frame() %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "expression"
  ) %>%
  mutate(
    group = factor(
      str_remove(sample, "-\\d+$"),
      levels = c("Normal","CIN1","CIN2","CIN3","Cancer")
    )
  )

qc_density <- ggplot(density_data, aes(x = expression, colour = group, group = sample)) +
  geom_density(alpha = 0.6, size = 0.5) +
  scale_colour_manual(values = group_colours) +
  labs(
    title = "QC Density Plot — GSE63514",
    subtitle = "Similar curve shapes = consistent normalisation",
    x = "Expression (log2)", y = "Density", colour = "Group"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "top"
  )

ggsave("figures/quality_control/QC_02_density.png", qc_density, width = 10, height = 6, dpi = 300)
print(qc_density)




# PCA PLOT

# Run PCA
pca_result <- prcomp(t(expr_matrix), scale. = TRUE, center = TRUE)

# How much variation does each PC explain?
pc_var <- round(summary(pca_result)$importance[2, 1:2] * 100, 1)
cat("PC1 explains:", pc_var[1], "% of variation\n")
cat("PC2 explains:", pc_var[2], "% of variation\n\n")

# Build plot data
pca_df <- data.frame(
  PC1   = pca_result$x[, 1],
  PC2   = pca_result$x[, 2],
  group = factor(
    str_remove(rownames(pca_result$x), "-\\d+$"),
    levels = c("Normal","CIN1","CIN2","CIN3","Cancer")
  )
)

qc_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = group, fill = group)) +
  stat_ellipse(alpha = 0.15, geom = "polygon", level = 0.8) +
  geom_point(size = 3, alpha = 0.8) +
  scale_colour_manual(values = group_colours) +
  scale_fill_manual(values = group_colours) +
  labs(
    title = "PCA Plot — GSE63514 Sample Clustering",
    subtitle = "Samples should cluster by disease grade",
    x = paste0("PC1 (", pc_var[1], "%)"),
    y = paste0("PC2 (", pc_var[2], "%)"),
    colour = "Group", fill = "Group"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "right"
  )

ggsave("figures/quality_control/QC_03_PCA.png", qc_pca, width = 10, height = 8, dpi = 300)
print(qc_pca)
