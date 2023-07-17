# Alpha diversity 
# SAB 7/17/2023

color_palette <- c("#0070FF", "#D75CE0", "#FF5EAA", "#FF8C76", "#FFC55A", "#F9F871")


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
varMF # p = 0.03, SIGNIFICANT

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

ggviolin(ps.meta, x = "OperatorGender", y = "Shannon$Shannon",
                  add = "boxplot", fill = "OperatorGender", palette = color_palette, title = "C", ylab = "Shannon's Diversity Index", xlab = "Gender") +
 theme(legend.position = "none")


