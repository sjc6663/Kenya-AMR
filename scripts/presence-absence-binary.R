# Presence/Absence Binary 
# SAB 1/29/24

# load packages
require(phyloseq)
require(microViz)
require(ggplot2)
require(vegan)
require(dplyr)
require(tibble)
require(viridis)
require(BiocManager)
require(microbiome)
require(ggpubr)
require(stringr)
require(patchwork)

# install.packages("remotes")
remotes::install_github("adrientaudiere/MiscMetabar")
library(MiscMetabar)

test <- as_binary_otu_table(ps, min_number = 1)

bin <- test

saveRDS(bin, file = "bovine-host-resistome/binary-ps.RDS")

# aggregate taxa ----
psbclass <- aggregate_taxa(bin, level = "Broadclass")

psbclass %>% plot_composition(average_by = "OperatorGender") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  labs(x = " ", y = "Relative Abundance", fill = "Resistance Type") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  ggtitle("A")

# merge samples by comparative categories (metadata) we want
mps <- merge_samples(bin, "OperatorGender")

# create a dataframe with the aggregated abundances for desired tax rank
phy <- mps %>% tax_glom(taxrank = "Broadclass") %>% psmelt()

# select only relevant columns
phy <- select(phy, c("Sample", "Abundance", "Broadclass"))

# get abundance as a percent and round to whole numbers
phy$percent <- phy$Abundance * 100
phy$percent <- round(phy$percent, digits = 0)

phy <- phy[order(phy$Broadclass), ]
phy
