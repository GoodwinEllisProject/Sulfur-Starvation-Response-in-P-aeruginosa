# Load the necessary package
library(dplyr)
library(readxl)
# Define the two datasets
dataset1 <- read.csv("Updated_table_LFC_1_RNAseq.csv", header = TRUE)
dataset2 <- read_xlsx("idmapping_2024_05_23.xlsx")

# Merge the datasets based on the ID variable
merged_rna_data <- full_join(dataset1, dataset2, by = c("GeneID" = "From"))
write.csv(merged_rna_data, "merged_rna_data.csv")

filtered_rna_dat <- merged_rna_data %>% 
  filter(abs(logFC) >= 2 & adj.P.Val <= 0.01)

# export as csv
write.csv(filtered_rna_dat, "rnaseq_strongDE.csv", row.names = FALSE)


# filter proteomics data
proteomics_dat <- read.csv("DE_results_proteomics.csv", header = TRUE)

filtered_proteo_dat <- proteomics_dat %>% 
  filter(abs(.[[3]]) >= 2 & .[[5]] <= 0.01)

# export as csv
write.csv(filtered_proteo_dat, "strongDE_prot.csv", row.names = FALSE)


# Assume df1 and df2 are your two dataframes with protein IDs in a column named "Protein_ID"

# Remove rows with missing values

filtered_proteo_dat <- na.omit(filtered_proteo_dat)
filtered_rna_dat <- na.omit(filtered_rna_dat)

# Find the common protein IDs
common_proteins <- intersect(filtered_rna_dat$Entry, filtered_proteo_dat$Protein.IDs)

# Create a pie chart
pie(table(c(filtered_rna_dat$Entry, filtered_proteo_dat$Protein.IDs)), main = "Protein ID Overlap")


