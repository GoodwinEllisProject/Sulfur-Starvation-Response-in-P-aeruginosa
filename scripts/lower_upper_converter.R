# Load necessary libraries
library(readr)

# Load data from CSV file
df <- read.csv("full_table_proteomics.csv")

# Convert first letter to uppercase
df$Gene.Name <- paste0(toupper(substr(df$Gene.Name, 1, 1)), substr(df$Gene.Name, 2, nchar(df$Gene.Name)))

# Save updated data to new CSV file
write_csv(df, "updated_fullde_proteomics.csv")
