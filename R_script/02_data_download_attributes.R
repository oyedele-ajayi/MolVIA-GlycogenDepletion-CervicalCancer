# DOWNLOAD GSE63514

# Load the download tool
library(GEOquery)

# Where to save the dataset
download_folder <- "data/raw"

# Download the dataset
gse63514 <- getGEO(
  GEO       = "GSE63514",      
  destdir   = download_folder,  
  GSEMatrix = TRUE,             
  AnnotGPL  = TRUE              
)

#  DATA EXPLORATION: FIRST LOOK AT THE DATASET

# we take the first item
gse_data <- gse63514[[1]]


# Samples (columns) 
cat("Number of samples:", ncol(gse_data), "\n")

# How many genes (rows)
cat("Number of genes/probes:", nrow(gse_data), "\n")

# What does the data structure look like?
cat("\nData class:", class(gse_data), "\n")

# Extract the sample information table
sample_info <- pData(gse_data)

# How many columns of information
cat("Number of information columns:", ncol(sample_info), "\n\n")

# Show the column names
print(colnames(sample_info))


# FIND THE DISEASE GROUP COLUMN

# These columns usually contain the group information in GEO
key_columns <- c("title", "source_name_ch1", 
                 "characteristics_ch1", 
                 "characteristics_ch1.1",
                 "characteristics_ch1.2",
                 "description")

# Show only columns that exist in the data
existing_cols <- key_columns[key_columns %in% colnames(sample_info)]

# Print the first 10 rows of these columns
print(sample_info[1:10, existing_cols])



# COUNT SAMPLES PER DISEASE GROUP

# which sample descriptions
print(sample_info$title[1:20])

print(table(sample_info$source_name_ch1))


# save the downloaded data object
save(gse63514, gse_data, sample_info,
     file = "data/raw/GSE63514_downloaded.RData")
