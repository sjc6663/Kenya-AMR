# Resistome - Beta Diversity 
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

ait <- psg %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "OperatorGender") %>% bdisp_get() # P = 0.35, ns

# test with PERMANOVA
mod3 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "HerdSize + HighestLevelEducation + OperatorGender + HerdSize*HighestLevelEducation + HighestLevelEducation*OperatorGender + HerdSize*HighestLevelEducation*OperatorGender",
    n_perms = 9999
  )

mod3 #P > 0.06, so ns, therefore we can remove the three-way interaction

mod2 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "HerdSize + HighestLevelEducation + OperatorGender + HerdSize*HighestLevelEducation + HighestLevelEducation*OperatorGender",
    n_perms = 9999
  )

mod2 #P > 0.10, so ns, therefore we can remove the two-way interactions

mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "HerdSize + HighestLevelEducation + OperatorGender",
    n_perms = 9999
  )

mod1 # Herd Size, R2 = 0.05, F(1,18) = 1.18, P = 0.11
      # Education, R2 = 0.14, F(3,18) = 1.07, P = 0.27
      # Gender, R2 = 0.05, F(1,18) = 1.13, P = 0.22, no significance for any group

C <- ait %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  #tax_transform(trans = "clr") %>%
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "OperatorGender", size = "HerdSize", axes = c(1,2), shape = "HighestLevelEducation") +
  scale_color_manual(values = short_color) +
  stat_ellipse(aes(group = OperatorGender, color = OperatorGender)) + 
  theme_classic() +
  labs(color = "Gender") +
  ggtitle("C") + 
 # labs(caption = "R2 = 0.04, F(1,22) = 0.96, P = 0.58, Padj = ___") +
  theme(text=element_text(size=15), #change font size of all text
        axis.text=element_text(size=20), #change font size of axis text
        axis.title=element_text(size=20), #change font size of axis titles
        plot.title=element_text(size=20), #change font size of plot title
        legend.text=element_text(size=20), #change font size of legend text
        legend.title=element_text(size=20)) #change font size of legend title 
C