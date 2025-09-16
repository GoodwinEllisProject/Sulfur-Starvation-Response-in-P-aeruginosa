# Load the necessary libraries
library(readr)
library(dplyr)

# Read the CSV file
data <- read_csv("full_table_proteomics.csv")

# Filter the data
filtered_data <- data %>%
  filter(data$significant == TRUE)

# Write the filtered data to a new CSV file
write_csv(filtered_data, "DE_prots.csv")

