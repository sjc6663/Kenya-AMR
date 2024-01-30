# Preliminary Kenya Analysis to Convert to R Markdown File - HERD SIZE (3 or less)
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
short_color <- c("#2e6f8f", "#83a561")


# Descriptive Stats --- Barplots ----

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# aggregate taxa ----
psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

#A <- psbclass %>% plot_composition(average_by = "HerdSize") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
#  scale_fill_manual(values = color_palette) + 
#  theme_bw() +
#  labs(x = " ", y = "Relative Abundance", fill = "Resistance Type") +
#  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
#  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n")) +
#  ggtitle("A")

#ggsave(A, file = "plots/size-relabund.tiff", dpi = 600)

# Descriptive Stats --- Relative Abundance Percentages ----

# merge samples by comparative categories (metadata) we want
mps <- merge_samples(psrel, "HerdSize")

# create a dataframe with the aggregated abundances for desired tax rank
phy <- mps %>% tax_glom(taxrank = "Broadclass") %>% transform_sample_counts(function(x) {x/sum(x)}) %>% psmelt()

# select only relevant columns
phy <- select(phy, c("Sample", "Abundance", "Broadclass"))

phy

## Plot with Percentages of Relative Abundance

phy$Sample <- factor(phy$Sample, levels = c("Less than 3", "More than 3"))


phy %>% 
  mutate(percent_labels = scales::percent(Abundance, accuracy = .01L)) %>% 
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
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  scale_x_discrete(labels = function(x) lapply(strwrap(x, width = 10, simplify = FALSE), paste, collapse="\n"))


ggsave(filename = "plots/size-relabund-percents.pdf", dpi = 600)

# Alpha Diversity --- Linear Regression and Violin Plots ----

adiv <- data.frame(
  "Shannon" = phyloseq::estimate_richness(ps, measures = "Shannon"),
  "MF" = phyloseq::sample_data(ps)$OperatorGender,
  "MilkSales" = phyloseq::sample_data(ps)$MilkSales,
  "Size" = phyloseq::sample_data(ps)$HerdSize,
  "Employees" = phyloseq::sample_data(ps)$NonFamilyEmployees,
  "LastABXUse" = phyloseq::sample_data(ps)$LastABXUse,
  "Route" = phyloseq::sample_data(ps)$ABXRoute,
  "Edu" = phyloseq::sample_data(ps)$HighestLevelEducation
)

mod <- lm(Shannon ~ Size, data = adiv)
summary(mod) # p = 0.03*

ps.meta <- meta(ps)
ps.meta$Shannon <- phyloseq::estimate_richness(ps, measures = "Shannon")

ps.meta$'' <- alpha(ps, index = 'shannon')

ggviolin(ps.meta, x = "size", y = "Shannon$Shannon",
         add = "boxplot", fill = "size", palette = short_color, title = " ", ylab = "Shannon's Diversity Index", xlab = " ") +
  theme(legend.position = "none") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  scale_y_continuous(limits = c(1, 7))

B <- plot_richness(ps, x="HerdSize", measures=c("Shannon"), color = "HerdSize") + 
  geom_boxplot() + 
  geom_jitter() +
  scale_color_manual(values = short_color) + 
  theme_classic() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(1,6.5)) +
  labs(x = " ", y = "Alpha Diversity Measure") +
  theme(text = element_text(size = 20), axis.title = element_text(size = 20)) +
  ggtitle("B")

#ggsave(B, file = "plots/size-alpha.tiff", dpi = 600)

# Beta Diversity --- PERMANOVA and PCoA ----

ait <- ps %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "HerdSize") %>% bdisp_get() # P = 0.64, ns

# test with PERMANOVA
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "HerdSize + HighestLevelEducation + OperatorGender + HerdSize*HighestLevelEducation + HerdSize*OperatorGender",
    n_perms = 9999
  )

mod1 #P > 0.05, so ns

mod2 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "HerdSize",
    n_perms = 9999
  )

mod2 # R2 = 0.04, F(1,22)  1.02, P = 0.32

ait %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  #tax_transform(trans = "clr") %>%
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "HerdSize", size = 6, axes = c(1,2)) +
  scale_color_manual(values = short_color) +
  stat_ellipse(aes(group = HerdSize, color = HerdSize)) + 
  theme_classic() +
  labs(color = "Herd Size") +
  ggtitle(" ") + 
  labs(caption = " ") +
  theme(text = element_text(size = 20)) 

# Differential Relative Abundance ----

psclass <- aggregate_taxa(ps, level = "Class")

# run AlDEx2 function
aldex2 <- ALDEx2::aldex(data.frame(phyloseq::otu_table(psclass)), phyloseq::sample_data(psclass)$HerdSize, taxa_rank = "all", norm = "CLR", method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 128, denom = "iqlr")

# look to see if anything is significant (this is for nonparametric, parametric use we.eBH)
sig_aldex2 <- aldex2 %>%
  filter(wi.eBH < 0.05)

sig_aldex2

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
