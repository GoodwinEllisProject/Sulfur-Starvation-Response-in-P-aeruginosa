#──────────────────────────────────────────────────────────────────────────────
# 0) INSTALL & LOAD PACKAGES
#──────────────────────────────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma","EnhancedVolcano"), ask=FALSE)

install.packages(c("tidyverse","ggpubr","plotly"), dependencies=TRUE)

library(limma)
library(EnhancedVolcano)
library(tidyverse)   # includes ggplot2, dplyr, tidyr
library(ggpubr)
library(plotly)

#──────────────────────────────────────────────────────────────────────────────
# 1) READ & PREPARE DATA
#──────────────────────────────────────────────────────────────────────────────
raw <- read.csv("Imputed_matrix.csv", header=TRUE, stringsAsFactors=FALSE)

# Capitalize Gene names
raw <- raw %>%
  mutate(Gene = paste0(toupper(substr(Gene,1,1)),
                       substr(Gene,2,nchar(Gene))))

# Keep a map for later
gene_map <- setNames(raw$Gene, raw$ProteinID)

# Build expression matrix
exprs <- raw[, 3:ncol(raw)] %>%
  set_names(sub(" .*", "", names(.)))  # drop suffix after space
rownames(exprs) <- raw$ProteinID

#──────────────────────────────────────────────────────────────────────────────
# 2) LIMMA DESIGN & CONTRASTS
#──────────────────────────────────────────────────────────────────────────────
# Define your 6 groups × 3 reps
sample_order <- rep(c("C2.5h","C5h","C7.5h","T2.5h","T5h","T7.5h"), each=3)
stopifnot(length(sample_order) == ncol(exprs))

group  <- factor(sample_order, 
                 levels = c("C2.5h","C5h","C7.5h","T2.5h","T5h","T7.5h"))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# All pairwise contrasts
contrasts <- makeContrasts(
  C2.5_vs_C5   = C2.5h - C5h,
  C2.5_vs_C7.5 = C2.5h - C7.5h,
  C2.5_vs_T2.5 = C2.5h - T2.5h,
  C2.5_vs_T5   = C2.5h - T5h,
  C2.5_vs_T7.5 = C2.5h - T7.5h,
  C5_vs_C7.5   = C5h   - C7.5h,
  C5_vs_T2.5   = C5h   - T2.5h,
  C5_vs_T5     = C5h   - T5h,
  C5_vs_T7.5   = C5h   - T7.5h,
  C7.5_vs_T2.5 = C7.5h - T2.5h,
  C7.5_vs_T5   = C7.5h - T5h,
  C7.5_vs_T7.5 = C7.5h - T7.5h,
  T2.5_vs_T5   = T2.5h - T5h,
  T2.5_vs_T7.5 = T2.5h - T7.5h,
  T5_vs_T7.5   = T5h   - T7.5h,
  levels = design
)

#──────────────────────────────────────────────────────────────────────────────
# 3) FIT & EXTRACT DE RESULTS
#──────────────────────────────────────────────────────────────────────────────
fit    <- lmFit(exprs, design)
fitC   <- contrasts.fit(fit, contrasts)
fitEB  <- eBayes(fitC)

DE.list <- lapply(colnames(contrasts), function(ct) {
  topTable(fitEB,
           coef        = ct,
           number      = Inf,
           adjust.method="BH",
           p.value     = 0.05,
           lfc         = 1) %>%
    tibble::rownames_to_column("ProteinID") %>%
    mutate(
      Gene     = gene_map[ProteinID],
      Contrast = ct
    )
})
DE.all <- bind_rows(DE.list)

# Save & summarize
write.csv(DE.all, "DE_proteins_all_contrasts.csv", row.names=FALSE)
message("DE table written: DE_proteins_all_contrasts.csv")
DE.all %>% count(Contrast, name="n_DE") %>% print()

#──────────────────────────────────────────────────────────────────────────────
# 4) VOLCANO PLOT FUNCTION (by contrast or by gene)
#──────────────────────────────────────────────────────────────────────────────
plot_volcano <- function(contrast_name) {
  df <- DE.all %>% filter(Contrast==contrast_name)
  EnhancedVolcano(df,
                  lab      = df$Gene,
                  x        = "logFC",
                  y        = "P.Value",
                  pCutoff  = 0.05,
                  FCcutoff = 1,
                  title    = contrast_name
  )
}

plot_volcano_by_gene <- function(gene_name, contrast_name) {
  df <- DE.all %>% filter(Gene==gene_name, Contrast==contrast_name)
  EnhancedVolcano(df,
                  lab      = df$Gene,
                  x        = "logFC",
                  y        = "P.Value",
                  pCutoff  = 0.05,
                  FCcutoff = 1,
                  title    = paste(gene_name, contrast_name)
  )
}

# Example:
plot_volcano("C2.5_vs_C5")
#plot_volcano_by_gene("BfrB","C2.5_vs_C5")

#──────────────────────────────────────────────────────────────────────────────
# 5) PIVOT TO LONG FORM & RECODE CONDITIONS
#──────────────────────────────────────────────────────────────────────────────
df_long2 <- raw %>%
  pivot_longer(
    cols          = ends_with("Intensity_imputed_intensity"),
    names_to      = "Sample",
    names_pattern = "^(.*)\\.Intensity_imputed_intensity$",
    values_to     = "Expression"
  ) %>%
  mutate(
    # original Sample like "control_2_5h_1"
    Condition = sub("_[0-9]+$", "", Sample),
    Condition = recode(Condition,
                       "control_2_5h" = "SC_2.5h",
                       "control_5h"   = "SC_5h",
                       "control_7_5h" = "SC_7.5h",
                       "test_2_5h"    = "SF_2.5h",
                       "test_5h"      = "SF_5h",
                       "test_7_5h"    = "SF_7.5h"
    ),
    # keep Gene capitalized
    Gene       = paste0(toupper(substr(Gene,1,1)),
                        substr(Gene,2,nchar(Gene)))
  )

#──────────────────────────────────────────────────────────────────────────────
# 6) PREPARE P‑VALUE ANNOTATIONS
#──────────────────────────────────────────────────────────────────────────────
#genes <- c("KatA","KatB","AhpC","AhpF","OhrR","Ohr", "SodA", "SodB", "PA3450")   # set your subset
#genes <- c("KatB","OhrR","Ohr", "SodA", "SodB", "PA3450")   # set your subset
#genes <- c("BfrB","PA0962", "PA4880")   # set your subset
#genes <- c("EdaA","Edd", "Pgl","Zwf")   # set your subset
#genes <- c("PA2603", "PA2602")   
genes <- c("Gor","PA4401", "PA3628","PA2821", "PA2813", "PA1890", "PA1655")   # set your subset


# compute expression ranges
expr_ranges <- df_long2 %>%
  filter(Gene %in% genes) %>%
  group_by(Gene) %>%
  summarize(
    maxExpr = max(Expression, na.rm=TRUE),
    range   = maxExpr - min(Expression, na.rm=TRUE)
  )

# build manual_results2
manual_results2 <- DE.all %>%
  filter(Gene   %in% genes,
         adj.P.Val <= 0.05) %>%
  separate(Contrast, into=c("g1","g2"), sep="_vs_") %>%
  mutate(
    group1 = recode(g1,
                    "C2.5"="SC_2.5h","C5"="SC_5h","C7.5"="SC_7.5h",
                    "T2.5"="SF_2.5h","T5"="SF_5h","T7.5"="SF_7.5h"
    ),
    group2 = recode(g2,
                    "C2.5"="SC_2.5h","C5"="SC_5h","C7.5"="SC_7.5h",
                    "T2.5"="SF_2.5h","T5"="SF_5h","T7.5"="SF_7.5h"
    ),
    p.adj.lab = case_when(
      adj.P.Val < 0.001 ~ "p < .001",
      TRUE              ~ paste0("p = ", sprintf("%.3f", adj.P.Val))
    )
  ) %>%
  left_join(expr_ranges, by="Gene") %>%
  group_by(Gene) %>%
  arrange(Gene, group1, group2) %>%
  mutate(
    step       = row_number() - 1,
    y.position = maxExpr + (0.05 + 0.10*step)*range
  ) %>%
  ungroup() %>%
  distinct(Gene, group1, group2, p.adj.lab, y.position)

#──────────────────────────────────────────────────────────────────────────────
# 7) DRAW THE FACETTED BOXPLOT + P‑VALUE BRACKETS
#──────────────────────────────────────────────────────────────────────────────
p_final <- df_long2 %>%
  filter(Gene %in% genes) %>%
  ggplot(aes(x = Condition, y = Expression, fill = Condition)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
  facet_wrap(~ Gene, scales = "free_y", ncol = 3) +   
  scale_fill_brewer(palette = "Set2") +
  labs(x = "Group", y = "log2(Protein Abundance)") +
  theme_bw(base_size = 14) +
  theme(
    legend.position  = "none",
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.title       = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  stat_pvalue_manual(
    data        = manual_results2,
    label       = "p.adj.lab",
    xmin        = "group1",
    xmax        = "group2",
    y.position  = "y.position",
    tip.length  = 0.02,
    size        = 3.0,
    inherit.aes = FALSE
  )

print(p_final)
ggsave("Glutathione_met.png", p_final, width = 9, height = 12)
