# Preliminary Kenya Analysis to Convert to R Markdown File - POST DIP USE
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

# load phyloseq object
ps <- readRDS("bovine-host-resistome/decontam-ps.RDS")

color_palette <- c("#81007F", "#D75CE0", "#FF5EAA", "#FF8C76", "#FFC55A", "#F9F871")
four_color <- c("#81007F", "#D75CE0", "#FF8C76", "#FFC55A")
short_color <- c("#D75CE0", "#FFC55A")

# Descriptive Stats --- Barplots ----

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# aggregate taxa ----
psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

psbclass %>% plot_composition(group_by = "PostDip") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = four_color)
# ggtitle("A")

# Descriptive Stats --- Relative Abundance Percentages ----

# merge samples by comparative categories (metadata) we want
mps <- merge_samples(psrel, "PostDip")

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
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "MF" = phyloseq::sample_data(ps)$OperatorGender,
  "MilkSales" = phyloseq::sample_data(ps)$MilkSales,
  "Size" = phyloseq::sample_data(ps)$HerdSize,
  "Employees" = phyloseq::sample_data(ps)$NonFamilyEmployees,
  "LastABXUse" = phyloseq::sample_data(ps)$LastABXUse,
  "Route" = phyloseq::sample_data(ps)$ABXRoute,
  "Edu" = phyloseq::sample_data(ps)$HighestLevelEducation,
  "dip" = phyloseq::sample_data(ps)$PostDip
)

mod <- lm(Shannon ~ dip, data = adiv)
summary(mod)

ps.meta <- meta(ps)
ps.meta$Shannon <- phyloseq::estimate_richness(ps, measures = "Shannon")

ps.meta$'' <- alpha(ps, index = 'shannon')

ggviolin(ps.meta, x = "PostDip", y = "Shannon$Shannon",
         add = "boxplot", fill = "PostDip", palette = short_color, title = " ", ylab = "Shannon's Diversity Index", xlab = " ") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(1, 7))

# Beta Diversity --- PERMANOVA and PCoA ----

ait <- ps %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "PostDip") %>% bdisp_get()

# test with PERMANOVA
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "PostDip",
    n_perms = 9999
  )

mod1 

ait %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  #tax_transform(trans = "clr") %>%
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "PostDip", size = 6, axes = c(1,2)) +
  scale_color_manual(values = short_color) +
  stat_ellipse(aes(group = PostDip, color = PostDip)) + 
  theme_classic() +
  labs(color = "Post Dip Use?") +
  ggtitle(" ") + 
  labs(caption = " ") +
  theme(text = element_text(size = 20)) 

# Differential Relative Abundance ----

# run AlDEx2 function
aldex2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(ps)), phyloseq::sample_data(ps)$PostDip, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
sig_aldex2 <- aldex2 %>%
  filter(wi.eBH < 0.05)

sig_aldex2

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
