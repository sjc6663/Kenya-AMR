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

set.seed(812)

# load phyloseq object
ps <- readRDS("data/decontam-ps-fixed.RDS")

color_palette <- c("#45337d", "#2e6f8f", "#218f8c", "#83a561", "#bfdb81", "#1fa187")
four_color <- c("#eae69e", "#bfdb81", "#83a561", "#48723e")
short_color <- c("#2e6f8f", "#83a561")

# Descriptive Stats --- Barplots ----

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")
# Descriptive Stats --- Relative Abundance Percentages ----

# merge samples by comparative categories (metadata) we want
mps <- merge_samples(psrel, "SampleBinary")

# create a dataframe with the aggregated abundances for desired tax rank
phy <- mps %>% tax_glom(taxrank = "Broadclass") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

# select only relevant columns
phy <- select(phy, c("Sample", "Abundance", "Broadclass"))

phy

## Plot with Percentages of Relative Abundance
phy %>% 
  mutate(percent_labels = scales::percent(Abundance)) %>% 
  ggplot(aes(x = Sample, y = Abundance, fill = Broadclass)) +
  geom_col() +
  geom_text(
    aes(label = percent_labels), 
    position = position_stack(vjust = 0.5),
    col = "white"
  ) + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  labs(x = " ", y = "Relative Abundance", fill = "Resistance Type") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20))

ggsave(filename = "plots/fixed/figure3-relabund-percents.jpg", dpi = 600)
