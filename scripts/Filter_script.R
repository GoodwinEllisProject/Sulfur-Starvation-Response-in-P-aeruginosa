# Read the CSV file into a data frame
df <- read.csv("Imputed_matrix.csv")

# Extract rows with "pvd" or "fpv" anywhere in the row
#pvd_fpv_rows <- df[grep("pvd|fpv", apply(df, 1, paste, collapse = ",")), ]
#pvd_rows <- df[grep("pvd", apply(df, 1, paste, collapse = ",")), ]
#biofilm <- df[grep("pil|fli|pil|las|rhl|las", apply(df, 1, paste, collapse = ",")), ]
phenazines_count <- df[grep("phz", apply(df, 1, paste, collapse = ",")), ]


# View the extracted rows
#print(pvd_rows)
#head (pvd_fpv_rows)
head(phenazines_count)

# save as csv
#write.csv(pvd_fpv_rows, "pvd_fpv.csv")
write.csv(phenazines_count, "phenazines_count.csv")

