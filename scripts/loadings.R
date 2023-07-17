# loadings for RDA plot
# SAB 07/17/2023

ps <- readRDS("data/decontam-ps.RDS")

met <- readxl::read_xlsx("kenya-metadata.xlsx",
                         sheet = "meta-numbers")

samp <- phyloseq::sample_data(met)
rownames(samp) <- met$`SampleID`

sample_data(ps) <- phyloseq(sample_data(samp))
sample_data(ps)

library(stats)
library(phyloseq)
library(microViz)
library(ggplot2)

psrel <- microbiome::transform(ps, "compositional")

sample_data(ps)$OperatorGender <- as.integer(sample_data(ps)$OperatorGender)
sample_data(ps)$MilkSales <- as.numeric(sample_data(ps)$MilkSales)
sample_data(ps)$HerdSize <- as.numeric(sample_data(ps)$HerdSize)
sample_data(ps)$NonFamilyEmployees <- as.numeric(sample_data(ps)$NonFamilyEmployees)
sample_data(ps)$PreDip <- as.numeric(sample_data(ps)$PreDip)
sample_data(ps)$PostDip <- as.numeric(sample_data(ps)$PostDip)
sample_data(ps)$LastABXUse <- as.numeric(sample_data(ps)$LastABXUse)
sample_data(ps)$ABXRoute <- as.numeric(sample_data(ps)$ABXRoute)
sample_data(ps)$HighestLevelEducation <- as.numeric(sample_data(ps)$HighestLevelEducation)


# alternatively, constrain variation on weight and female
constrained_aitchison_rda <- psrel %>%
  tax_transform("clr") %>%
  ord_calc(constraints = c("MilkSales", "HerdSize", "NonFamilyEmployees", "PreDip", "PostDip", "LastABXUse", "ABXRoute", "HighestLevelEducation"), 
           method = "RDA") # constraints --> RDA

constrained_aitchison_rda %>%
  ord_plot(colour = "OperatorGender", constraint_vec_length = 2)
