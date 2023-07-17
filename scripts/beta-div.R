# beta diversity 
# SAB 07/17/2023

# load packages
library(phyloseq)
library(microViz)
library(vegan)
library(ggplot2)

# color scheme: Viridis mako / microshades micro_cvd_blue
color_palette <- c("#0070FF", "#D75CE0", "#FF5EAA", "#FF8C76", "#FFC55A", "#F9F871")

# load phyloseq object
ps <- readRDS("data/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")


# clr transform phyloseq objects
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()


dist_mat <- phyloseq::distance(transps, method = "euclidean")

vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$OperatorGender) # p = 0.845, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$HerdSize) # p = 0.086, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$NonFamilyEmployees) # p = 0.88, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$LastABXUse) # p = 0.649, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$ABXRoute) # p = 0.928, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$HighestLevelEducation) # p = 0.54, ns

# PCA Plots ------------------------------------------------------------------------------------------------------

psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "OperatorGender") +
  scale_color_manual(values = color_palette) +
  # stat_ellipse(aes(group = OperatorGender, color = OperatorGender)) + 
  theme_classic() +
  ggtitle("A") + 
  labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")
