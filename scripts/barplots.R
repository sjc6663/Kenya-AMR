# Barplots for Descriptive Stats
# SAB 07/17/2023

# color scheme: Viridis mako / microshades micro_cvd_blue
color_palette <- c("#81007F", "#D75CE0", "#FF5EAA", "#FF8C76", "#FFC55A", "#F9F871")
short_color <- c("#D75CE0", "#FFC55A")

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)
library(BiocManager)
library(stringr)
library(microshades)
library(microbiome)

# read in phyloseq object
ps <- readRDS("bovine-host-resistome/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# aggregate taxa ----
psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

A <- psbclass %>% plot_composition(group_by = "OperatorGender") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("A")

sample_data(psbclass)$HighestLevelEducation <- factor(sample_data(psbclass)$HighestLevelEducation,
                                                            levels = c("Primary Level", "Secondary Level", "Tertiary eduction", "Diploma/University degree"))


B <- psbclass %>% plot_composition(group_by = "HighestLevelEducation", sample.sort = "OperatorGender", x.label = "OperatorGender") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("B")

B

ggsave(filename = "plots/relabund-gender-edu.pdf", dpi = 600, width = 10, height = 8)
