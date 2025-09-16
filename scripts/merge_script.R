
# This is a short merg script

library(tidyverse)
# Read the CSV file into a data frame
df1 <- read.csv("phenazines.csv")
df2 <- read.csv("phenazines_count.csv")


merged_df <- left_join(df2,df1, by = c("ProteinID" = "Gene.Name"))

#write.csv(pvd_fpv_rows, "pvd_fpv.csv")
write.csv(merged_df, "DE_count_phz_prot.csv")

