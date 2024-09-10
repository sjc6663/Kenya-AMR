# Resistome - Alpha Diversity 
# SAB 2/24/24

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

# viridis, palette option "viridis"

color_palette <- c("#45337d", "#2e6f8f", "#218f8c", "#83a561", "#bfdb81", "#1fa187")
short_color <- c("#45337d", "#83a561")

# glom everything at the gene level so we don't have duplicates
psg <- tax_glom(ps, taxrank = "Gene")

adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(psg, measures = "Shannon"),
  "G" = phyloseq::sample_data(psg)$OperatorGender,
  "HS" = phyloseq::sample_data(psg)$HerdSize,
  "E" = phyloseq::sample_data(psg)$HighestLevelEducation
)

adiv

# linear model, we know that education can have an impact on herd size, so we need to look at that interaction
# we also know that gender can have an impact on education level, so we need to look at that interaction
# we don't anticipate an interaction between herd size and gender, however, we could possibly see a third level interaction between all three groups
mod3 <- lm(Shannon ~ G + E + HS + G*E + HS*E + G*E*HS, data = adiv)

summary(mod3) # the p-value of our three-way interaction is 0.75, so ns == we can remove it from the model

mod2 <- lm(Shannon ~ G + E + HS + G*E + HS*E, data = adiv)

summary(mod2) # the p-value of our 2 way interactions are all greater than 0.05 so ns == we can remove them from the model

mod <- lm(Shannon ~ G + E + HS, data = adiv)

summary(mod) # the p-values of our variables of interest are all greater than 0.05, so ns
# Gender, P = 0.91
# Education, P = , this however may need an anova to check for sure since it's more than 2 variables
# Herd Size, P = 0.21

## Boxplot ----
B <- plot_richness(ps, x="HerdSize", measures=c("Shannon"), shape = "OperatorGender", color = "HighestLevelEducation") + 
  geom_boxplot() + 
  geom_jitter(size = 3) +
  scale_color_manual(values = color_palette) + 
  theme_classic() +
  #theme(legend.position = "none") +
  scale_y_continuous(limits = c(1,6)) +
  labs(x = "Herd Size", y = "Shannon's Diversity Index", color = "Gender", shape = "Education") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  ggtitle("A") +
  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n"))
B
