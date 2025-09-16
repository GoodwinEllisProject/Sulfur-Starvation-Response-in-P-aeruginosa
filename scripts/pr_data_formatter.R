# Load required libraries
library(tidyverse)

process_plate_data <- function(file_path) {
  # Read the data
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Find all rows containing "Value"
  value_rows <- which(data[,1] == "Value")
  
  # Initialize list to store results from each reading
  all_readings <- list()
  
  # Process each reading
  for(i in seq_along(value_rows)) {
    # Extract 4 rows after each "Value" row (A-F)
    start_row <- value_rows[i] + 1
    end_row <- start_row + 3
    
    # Extract values for this reading
    reading_values <- data[start_row:end_row, 2:6]  # Columns 2-6 contain the numeric values
    
    # Convert to numeric, removing any non-numeric characters
    reading_values <- apply(reading_values, 2, function(x) as.numeric(gsub("[^0-9.]", "", x)))
    
    # Add row names
    rownames(reading_values) <- c("Control", "Control_Dye", "test", "test_Dye")
    
    # Store in list
    all_readings[[i]] <- reading_values
  }
  
  # Calculate summary statistics for each reading
  summary_stats <- lapply(seq_along(all_readings), function(i) {
    reading <- all_readings[[i]]
    data.frame(
      Reading = i,
      Condition = rownames(reading),
      Mean = rowMeans(reading, na.rm = TRUE),
      SD = apply(reading, 1, sd, na.rm = TRUE)
    )
  })
  
  # Combine all summaries
  final_summary <- do.call(rbind, summary_stats)
  
  # Save results
  write.csv(final_summary, "analyzed_results.csv", row.names = FALSE)
  
  # Print summary
  print("Summary of all readings:")
  print(final_summary)
  
  return(list(raw_values = all_readings, summary = final_summary))
}

# Use the function
file_path <- "Labile_Fe_2.5h.csv"  # Replace with your file path
results <- process_plate_data(file_path)

# Load the data
data <- read.csv("analyzed_results.csv")

# Convert the data to long format
library(tidyr)
data_long <- pivot_longer(data, cols = c(Mean, SD), names_to = "Metric")

# Create the line plot
library(ggplot2)
ggplot(data_long, aes(x = Reading, y = value, color = Condition)) + 
  geom_line() + 
  facet_wrap(~ Metric, scales = "free_y") + 
  theme_classic()

# Export the mean values in CSV format
write.csv(data[, c("Condition", "Mean")], "mean_values.csv", row.names=FALSE)

# Load required libraries
library(tidyr)
library(dplyr)
library(readr)

# Read the CSV file
data <- read_csv("analyzed_results.csv")

# Reshape data to wide format for GraphPad Prism
prism_data <- data %>%
  select(Reading, Condition, Mean) %>%        # Focus on relevant columns
  pivot_wider(names_from = Condition,         # Pivot conditions into columns
              values_from = Mean,             # Use Mean values
              values_fill = NA)               # Fill missing values with NA

# Write the reformatted table to a new CSV
write_csv(prism_data, "formatted_for_prism.csv")

# Preview the output
print(head(prism_data))



