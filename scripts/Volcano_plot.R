#Adapted from https://sdgamboa.github.io/post/2020_volcano/

# load necessary packages
library(limma)
library(ggplot2)
library(tidyverse)
library(ggrepel)
library(kableExtra)
library(ggpubr)

# read in gene expression data
gene_exp <- read.csv("full DE_table.csv")

# replace missing gene names with corresponding gene ID
gene_exp$gene.name <- ifelse(is.na(gene_exp$gene.name), gene_exp$GeneID, gene_exp$gene.name)

# A short function for outputting the tables
knitr_table <- function(x) {
  x %>% 
    knitr::kable(format = "html", digits = Inf, 
                 format.args = list(big.mark = ",")) %>%
    kableExtra::kable_styling(font_size = 15)
}

head(gene_exp) %>% 
  knitr_table()


p1 <- ggplot(gene_exp, aes(logFC, -log(adj.P.Val,10))) + # -log10 conversion  
  geom_point(size = 2/5) +
  xlab(expression("log"[2]*"FC")) + 
  ylab(expression("-log"[10]*"adj.P.val"))
p1

# adding color to differntially expressed genes
data <- gene_exp %>% 
  mutate(
    Expression = case_when(logFC >= 1 &  adj.P.Val <= 0.05 ~ "Up-regulated",
                           logFC <= -1 & adj.P.Val <= 0.05 ~ "Down-regulated",
                           TRUE ~ "Unchanged")
  )
head(data) %>% 
  knitr_table()

#count genes
data %>% 
  count(Expression) %>% 
  knitr_table()

# Add colors 
# Customize the axes to be bold, in Arial, with a font size of 16 for Dr. Ellis
p2 <- ggplot(data, aes(logFC, -log(adj.P.Val, 10))) +
  theme_pubr() + 
  geom_point(aes(color = Expression), size = 2) +
  xlab(expression(bold("LogFC"))) +  # Make the x-axis label bold
  ylab(expression(bold(-log[10] ~ "(FDR adjusted p-value)"))) +  # Make the y-axis label bold
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
  guides(colour = "none") +  # Remove the expression legend
  theme(
    axis.text = element_text(face = "bold", family = "Arial", size = 16),  # Make axis text bold
    axis.title = element_text(face = "bold", family = "Arial", size = 18)  # Make axis labels bold
  )
p2



# p2 <- ggplot(data, aes(logFC, -log(adj.P.Val,10))) +
#   theme_pubr() + 
#   geom_point(aes(color = Expression), size = 2) +
#   xlab(expression("LogFC")) + 
#   ylab(expression(-log[10] ~ "(FDR adjusted p-value)")) +
#   scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) +
#   guides(colour = guide_legend(override.aes = list(size=1.5))) 
# p2

# Since we already know that the genes towards the right are up-regulated and the genes 
# towards the left are down-regulated, it would be more informative if we colored the points 
# according to their significance level instead. Let’s create another column, named ‘Significance’, and 
# classify the genes according to significance thresholds (0.05, 0.01, and 0.001):

data <- data %>% 
  mutate(
    Significance = case_when(
      abs(logFC) >= 1 & adj.P.Val <= 0.05 & adj.P.Val > 0.01 ~ "adj.P.Val 0.05", 
      abs(logFC) >= 1 & adj.P.Val <= 0.01 & adj.P.Val > 0.001 ~ "adj.P.Val 0.01",
      abs(logFC) >= 1 & adj.P.Val <= 0.001 ~ "adj.P.Val 0.001", 
      TRUE ~ "Unchanged")
  )
head(data) %>% 
  knitr_table()

p3 <- ggplot(data, aes(logFC, -log(adj.P.Val,10))) +
  theme_pubr() + 
  geom_point(aes(color = Significance), size = 1) +
  xlab(expression("logFC")) + 
  ylab(expression("-log"[10]*"adj.P.Val")) +
  scale_color_viridis_d() +
  guides(colour = guide_legend(override.aes = list(size=1.5))) 

p3
ggsave("funny.png", width = 7, height = 8, dpi = 400)

data %>% 
  count(Expression, Significance) %>% 
  knitr_table()

#adding labels

top <- 10
top_genes <- bind_rows(
  data %>% 
    filter(Expression == 'Up-regulated') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(top),
  data %>% 
    filter(Expression == 'Down-regulated') %>% 
    arrange(adj.P.Val, desc(abs(logFC))) %>% 
    head(top)
)
top_genes %>% 
  knitr_table()

p4 <-  p2 +
  geom_label_repel(data = top_genes,
                   mapping = aes(logFC, -log(adj.P.Val,10), label = gene.name),
                   size = 3)
p4
ggsave("Volcano_plot_2.png", width = 8.9, height = 5, dpi = 600)
