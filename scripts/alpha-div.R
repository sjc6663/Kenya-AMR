# Alpha diversity 
# SAB 7/17/2023

color_palette <- c("#81007F", "#D75CE0", "#FF5EAA", "#FF8C76", "#FFC55A", "#F9F871")
short_color <- c("#D75CE0", "#FFC55A")

set.seed(81299)

library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)
library(BiocManager)
library(stringr)
library(microbiome)
library(ggpubr)
library(knitr)
library(dplyr)
library(car)
library(misty)
library(patchwork)

# read in phyloseq object
ps <- readRDS("data/decontam-ps.RDS")

# create data frame ----
# create data frame with relevant metadata for comparison
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

# test variances ---------------------------------------------------------------------------------------------------
varMF <- var.test(Shannon ~ MF, data = adiv, 
                  alternative = "two.sided")
varMF # p = 0.03*, SIGNIFICANT

varMS <- leveneTest(Shannon ~ MilkSales, data = adiv) # have to do levene's test because its more than 2 groups
varMS # p = 0.53, not sig

varS <- leveneTest(Shannon ~ Size, data = adiv)
varS # p = 0.21, not sig

varE <- var.test(Shannon ~ Employees, data = adiv, 
                  alternative = "two.sided")
varE # p = 0.03, SIGNIFICANT

varABX <- leveneTest(Shannon ~ LastABXUse, data = adiv)
varABX # p = 0.59, not sig

varR <- leveneTest(Shannon ~ Route, data = adiv)
varR # p = 0.75, not sig

varEDU <- leveneTest(Shannon ~ Edu, data = adiv)
varEDU # p = 0.54, not sig

# since we have some unequal variances we need to run Welch's t-test for our statistical analysis

# Welch's T-test --------------------------------------------------------------------------------------------------

# male female
wtestMF <- oneway.test(Shannon ~ MF, data = adiv, var.equal = FALSE)
wtestMF # p = 0.58 ns

# herd size
wtestMS <- oneway.test(Shannon ~ Size, data = adiv, var.equal = FALSE)
wtestMS # p = 0.20, ns

# employees
wtestE <- oneway.test(Shannon ~ Employees, data = adiv, var.equal = FALSE)
wtestE # p = 0.70, ns

# last ABX use
wtestABX <- oneway.test(Shannon ~ LastABXUse, data = adiv, var.equal = FALSE)
wtestABX # p = 0.74, ns

# eduation
wtestEDU <- oneway.test(Shannon ~ Edu, data = adiv, var.equal = FALSE)
wtestEDU # p = 0.95, ns

# Violin Plots -----------------------------------------------------------------------------------------------
# https://microbiome.github.io/tutorials/PlotDiversity.html 
ps.meta <- meta(ps)
ps.meta$Shannon <- phyloseq::estimate_richness(ps, measures = "Shannon")

ps.meta$'' <- alpha(ps, index = 'shannon')

# Male Female ------------------------------------

v_A <- ggviolin(ps.meta, x = "OperatorGender", y = "Shannon$Shannon",
                add = "boxplot", fill = "OperatorGender", palette = c("#D75CE0", "#FFC55A"), title = "A", ylab = "Shannon's Diversity Index", xlab = "Gender")  +
  theme(legend.position = "none")
v_A

# create a list of pairwise comaprisons
MF <- unique(adiv$MF) # get the variables

MF_pair <- combn(seq_along(MF), 2, simplify = FALSE, FUN = function(i)MF[i])


v_A <- v_A + stat_compare_means(comparisons = MF_pair, label = "p.signif",
                                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_A

# Non-Family Employees

v_B <- ggviolin(ps.meta, x = "NonFamilyEmployees", y = "Shannon$Shannon",
                add = "boxplot", fill = "NonFamilyEmployees", palette = c("#D75CE0", "#FFC55A"), title = "B", ylab = "Shannon's Diversity Index", xlab = "Non-Family Employees")  +
  theme(legend.position = "none")
v_B

# create a list of pairwise comaprisons
NFE <- unique(adiv$Employees) # get the variables

NFE_pair <- combn(seq_along(NFE), 2, simplify = FALSE, FUN = function(i)NFE[i])


v_B <- v_B + stat_compare_means(comparisons = NFE_pair, label = "p.signif",
                                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_B

# Herd Size
v_C <- ggviolin(ps.meta, x = "HerdSize", y = "Shannon$Shannon",
                add = "boxplot", fill = "HerdSize", palette = c("#D75CE0", "#FFC55A", "#F9F871"), title = "C", ylab = "Shannon's Diversity Index", xlab = "Herd Size")  +
  theme(legend.position = "none")
v_C

# create a list of pairwise comaprisons
HS <- unique(adiv$Size) # get the variables

HS_pair <- combn(seq_along(HS), 2, simplify = FALSE, FUN = function(i)HS[i])


v_C <- v_C + stat_compare_means(comparisons = HS_pair, label = "p.signif",
                                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_C

# Last ABX Use
v_D <- ggviolin(ps.meta, x = "LastABXUse", y = "Shannon$Shannon",
                add = "boxplot", fill = "LastABXUse", palette = c("#D75CE0", "#FFC55A", "#FF8C76", "#F9F871"), title = "D", ylab = "Shannon's Diversity Index", xlab = "Last Antibiotic Use")  +
  scale_x_discrete(limits = c("Last year", "Last month", "This month", "This week")) +
   theme(legend.position = "none")
v_D

# create a list of pairwise comaprisons
ABX <- unique(adiv$LastABXUse) # get the variables

ABX_pair <- combn(seq_along(ABX), 2, simplify = FALSE, FUN = function(i)ABX[i])


v_D <- v_D + stat_compare_means(comparisons = ABX_pair, label = "p.signif",
                                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_D

# Highest Edu Level

v_E <- ggviolin(ps.meta, x = "HighestLevelEducation", y = "Shannon$Shannon",
                add = "boxplot", fill = "HighestLevelEducation", palette = c("#D75CE0", "#FFC55A", "#FF8C76", "#F9F871"), title = "E", ylab = "Shannon's Diversity Index", xlab = "Highest Education Level")  +
  scale_x_discrete(limits = c("Primary Level", "Secondary Level", "Tertiary eduction", "Diploma/University degree")) +
  theme(legend.position = "none")
v_E

# create a list of pairwise comaprisons
edu <- unique(adiv$Edu) # get the variables

edu_pair <- combn(seq_along(edu), 2, simplify = FALSE, FUN = function(i)edu[i])


v_E <- v_E + stat_compare_means(comparisons = edu_pair, label = "p.signif",
                                symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
)
v_E

(v_A/v_B/v_C)|(v_D/v_E)

ggsave(filename = "plots/alpha-all.pdf", dpi = 600, width = 15, height = 20)
