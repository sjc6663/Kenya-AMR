# Preliminary Kenya Analysis to Convert to R Markdown File - GENDER
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
ps <- readRDS("data/decontam-ps.RDS")

color_palette <- c("#45337d", "#2e6f8f", "#218f8c", "#83a561", "#bfdb81", "#1fa187")
four_color <- c("#eae69e", "#bfdb81", "#83a561", "#48723e")
short_color <- c("#2e6f8f", "#83a561") "#2db27d",

# Descriptive Stats --- Barplots ----

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# aggregate taxa ----
psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

A <- psbclass %>% plot_composition(average_by = "OperatorGender") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  theme_bw() +
  labs(x = " ", y = "Relative Abundance", fill = "Resistance Type") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  ggtitle("A")
A

#ggsave(A, file = "bovine-host-resistome/plots/gender-relabund.pdf", dpi = 600)

# Descriptive Stats --- Relative Abundance Percentages ----

# merge samples by comparative categories (metadata) we want
mps <- merge_samples(psrel, "OperatorGender")

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
  "Shannon" = phyloseq::estimate_richness(ps, measures = c("Shannon", "Observed")),
  "MF" = phyloseq::sample_data(ps)$OperatorGender,
  "MilkSales" = phyloseq::sample_data(ps)$MilkSales,
  "Size" = phyloseq::sample_data(ps)$HerdSize,
  "Employees" = phyloseq::sample_data(ps)$NonFamilyEmployees,
  "LastABXUse" = phyloseq::sample_data(ps)$LastABXUse,
  "Route" = phyloseq::sample_data(ps)$ABXRoute,
  "Edu" = phyloseq::sample_data(ps)$HighestLevelEducation
)

mod <- lm(Shannon.Observed ~ MF + Size + MF*Size, data = adiv)
summary(mod) # P-value = 0.43, ns (interaction is 0.49 so also not significant)

ps.meta <- meta(ps)
ps.meta$Shannon <- phyloseq::estimate_richness(ps, measures = "Shannon")

ps.meta$'' <- alpha(ps, index = 'shannon')

ggviolin(ps.meta, x = "OperatorGender", y = "Shannon$Shannon",
              add = "boxplot", fill = "OperatorGender", palette = short_color, title = " ", ylab = "Shannon's Diversity Index", xlab = " ") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(1, 7))

B <- plot_richness(ps, x="OperatorGender", measures=c("Observed"), color = "OperatorGender") + 
  geom_boxplot() + 
  geom_jitter() +
  scale_color_manual(values = short_color) + 
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(1,800)) +
  labs(x = " ", y = "Alpha Diversity Measure") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  ggtitle("B")

ggsave(file = "bovine-host-resistome/plots/gender-alpha-observed.tiff", dpi = 600)

# Beta Diversity --- PERMANOVA and PCoA ----

ait <- ps %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "OperatorGender") %>% bdisp_get() # P = 0.32

# test with PERMANOVA
# first check if there is an interaction, if none proceed
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "OperatorGender + HerdSize + OperatorGender*HerdSize",
    n_perms = 9999
  )

mod1 #P = 0.33, no interaction

mod2 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "OperatorGender",
    n_perms = 9999
  )

mod2 # R2 = 0.04, F(1,22) = 0.96, P = 0.58, Padj = ___

C <- ait %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  #tax_transform(trans = "clr") %>%
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "OperatorGender", size = 6, axes = c(1,2)) +
  scale_color_manual(values = short_color) +
  stat_ellipse(aes(group = OperatorGender, color = OperatorGender)) + 
  theme_classic() +
  labs(color = "Gender") +
  ggtitle("C") + 
  labs(caption = "R2 = 0.04, F(1,22) = 0.96, P = 0.58, Padj = ___") +
  theme(text=element_text(size=15), #change font size of all text
          axis.text=element_text(size=20), #change font size of axis text
          axis.title=element_text(size=20), #change font size of axis titles
          plot.title=element_text(size=20), #change font size of plot title
          legend.text=element_text(size=20), #change font size of legend text
          legend.title=element_text(size=20)) #change font size of legend title   

#ggsave(C, file = "bovine-host-resistome/plots/gender-beta.pdf", dpi = 600)

# Differential Relative Abundance ----

# run AlDEx2 function
psbclass <- aggregate_taxa(ps, level = "Class")

aldex2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(psbclass)), phyloseq::sample_data(psbclass)$OperatorGender, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
sig_aldex2 <- aldex2 %>%
  filter(wi.eBH < 0.05)

sig_aldex2

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)

all2 <- ggarrange(A,B,C, ncol = 1, nrow = 3)
all2

ggsave(all2, file = "bovine-host-resistome/plots/gender-combined.pdf", dpi = 600, width = 6, height = 12)
