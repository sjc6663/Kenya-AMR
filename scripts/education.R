# Preliminary Kenya Analysis to Convert to R Markdown File - EDUCATION
# SAB 12/5/2023

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
ps <- readRDS("bovine-host-resistome/decontam-ps.RDS")

color_palette <- c("#45337d", "#2e6f8f", "#218f8c", "#83a561", "#bfdb81", "#1fa187")
four_color <- c("#eae69e", "#bfdb81", "#83a561", "#48723e")
short_color <- c("#2e6f8f", "#83a561")

# Descriptive Stats --- Barplots ----

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# aggregate taxa ----
psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

sample_data(psbclass)$HighestLevelEducation <- factor(sample_data(psbclass)$HighestLevelEducation,
                                                      levels = c("Primary Level", "Secondary Level", "Tertiary education", "Diploma or University degree"))


A <- psbclass %>% plot_composition(average_by = "HighestLevelEducation") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) + 
  theme_bw() +
  labs(x = " ", y = "Relative Abundance", fill = "Resistance Type") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) +
  ggtitle("A")
A

#ggsave(A, file = "bovine-host-resistome/plots/education-relabund.pdf", dpi = 600)

# Descriptive Stats --- Relative Abundance Percentages ----

# merge samples by comparative categories (metadata) we want
mps <- merge_samples(psrel, "HighestLevelEducation")

# create a dataframe with the aggregated abundances for desired tax rank
phy <- mps %>% tax_glom(taxrank = "Broadclass") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

# select only relevant columns
phy <- select(phy, c("Sample", "Abundance", "Broadclass"))

# get abundance as a percent and round to whole numbers
phy$percent <- phy$Abundance * 100
phy$percent <- round(phy$percent, digits = 0)

phy

# Alpha Diversity --- Linear Regression and Violin Plots ----

adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = c("Shannon")),
  "MF" = phyloseq::sample_data(ps)$OperatorGender,
  "MilkSales" = phyloseq::sample_data(ps)$MilkSales,
  "Size" = phyloseq::sample_data(ps)$HerdSize,
  "Employees" = phyloseq::sample_data(ps)$NonFamilyEmployees,
  "LastABXUse" = phyloseq::sample_data(ps)$LastABXUse,
  "Route" = phyloseq::sample_data(ps)$ABXRoute,
  "Edu" = phyloseq::sample_data(ps)$HighestLevelEducation
)

# test normality to ensure we can use anova
shapiro.test(adiv$Shannon.Shannon) # P = 0.0026, NOT NORMAL DATA

# data is not normal so we need to do kruskal wallis instead of ANOVA
kruskal.test(Shannon.Shannon ~ Edu, data = adiv) # P = 0.78, ns

sample_data(ps)$HighestLevelEducation <- factor(sample_data(ps)$HighestLevelEducation,
                                                      levels = c("Primary Level", "Secondary Level", "Tertiary education", "Diploma or University degree"))


B <- plot_richness(ps, x="HighestLevelEducation", measures=c("Shannon"), color = "HighestLevelEducation") + 
  geom_boxplot() + 
  geom_jitter() +
  scale_color_manual(values = color_palette) + 
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(1,6)) +
  labs(x = " ", y = "Shannon's Diversity Index") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  ggtitle("B") +
  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n"))

#ggsave(B, file = "bovine-host-resistome/plots/education-alpha.tiff", dpi = 600)

# Beta Diversity --- PERMANOVA and PCoA ----

ait <- ps %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "HighestLevelEducation") %>% bdisp_get() # P > 0.1 for all values

# test with PERMANOVA for confounding variable interaction
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "HighestLevelEducation + OperatorGender + HerdSize + HighestLevelEducation*OperatorGender + HighestLevelEducation*HerdSize",
    n_perms = 9999
  )

mod1 # P > 0.1 so there is no effect due to interaction

mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "HighestLevelEducation",
    n_perms = 9999
  )

mod1 # R2 = 0.12, F(3, 20) = 0.93, P = 0.66, ns

ait %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  #tax_transform(trans = "clr") %>%
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "HighestLevelEducation", size = 6, axes = c(1,2)) +
  scale_color_manual(values = four_color) +
  stat_ellipse(aes(group = HighestLevelEducation, color = HighestLevelEducation)) + 
  theme_classic() +
  labs(color = "Education") +
  ggtitle(" ") + 
  labs(caption = " ") +
  theme(text = element_text(size = 20)) 
