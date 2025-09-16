# Read the CSV file into a data frame
df <- read.csv("count_table.csv")

# Extract rows with "pvd" or "fpv" anywhere in the row
#pvd_fpv_rows <- df[grep("pvd|fpv", apply(df, 1, paste, collapse = ",")), ]
#pvd_rows <- df[grep("pvd", apply(df, 1, paste, collapse = ",")), ]
#biofilm <- df[grep("pil|fli|pil|las|rhl|las", apply(df, 1, paste, collapse = ",")), ]
oxidative_stress <- df[grep("lsf|ohr|sodB", apply(df, 1, paste, collapse = ",")), ]


# View the extracted rows
#print(pvd_rows)
#head (pvd_fpv_rows)
head(oxidative_stress)

# save as csv
#write.csv(pvd_fpv_rows, "pvd_fpv.csv")
write.csv(oxidative_stress, "oxidativestress_count.csv")

