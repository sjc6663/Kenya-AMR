# Decontam
## SAB 11/16/2023

# load packages ----
library(BiocManager)
library(stringr)
library(tidyverse)
library(tidylog)
library(phyloseq)
library(vegan) #version 2.6-4
library(ggplot2)
library(decontam)
library(readxl)
library(microViz)

# load ps objects and other data ----
countsDF <- read.delim("data/countmatrix-cleanedall.txt", sep = "\t") 
met <- read_excel("data/kenya-metadata.xlsx")
genes <- read.delim("data/geneinfo-all.txt", sep = "\t") %>%
  select(-c(MEG_ID)) %>%
  unique()

ps <- readRDS("data/rawps.rds")

# dummy code samples vs controls ----
ps <- ps %>% 
  ps_mutate(
    SampleBinary = if_else(str_detect(OperatorGender,"Control"), true = "Control", false = "Sample")
  ) 

unknowns <- c("24-3", "2-3S", "23-1", "24-1")

ps <- ps %>% 
  subset_samples(
    FarmID != "24-3"
  )

ps <- ps %>% 
  subset_samples(
    FarmID != "2-3S"
  )

ps <- ps %>% 
  subset_samples(
    FarmID != "23-1"
  )

ps <- ps %>% 
  subset_samples(
    FarmID != "24-1"
  )


ps <- ps %>% 
  ps_mutate(SampleBinary = case_when(
    SampleBinary == "Control" & str_detect(OperatorGender, "Pos") ~ "Sample",
    SampleBinary == "Control" & !str_detect(OperatorGender, "Pos") ~ "Control",
    SampleBinary == "Sample" ~ "Sample"
  ))


# inspect # of counts per sample
sam <- as.data.frame(sample_data(ps))
sam$GeneCounts <- sample_sums(ps)
sam <- sam[order(sam$GeneCounts), ]
sam$Index <- seq(nrow(sam))


# plot to look at it 
ggplot(sam, aes(x = Index, y = GeneCounts, color = SampleBinary)) +
  geom_point()

## Negative Control ----
sample_data(ps)$is.neg <- sample_data(ps)$SampleBinary == "Control"
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)

# make a phyloseq object of only the contaminants and controls
con <- rownames(contamdf.prev[contamdf.prev$contaminant == "TRUE",])
conps <- prune_taxa(rownames(ps@tax_table) %in% con, ps)
negps <- subset_samples(conps, is.neg == "TRUE")

# remove contaminant sequences
contams <- rownames(contamdf.prev[contamdf.prev$contaminant == "TRUE",])
contams #135 gene contaminants

# remove contaminant sequences
nocontam <- prune_taxa(!rownames(ps@tax_table) %in% contams, ps)

saveRDS(nocontam, "bovine-host-resistome/decontam-ps.rds")

## Positive Control ----

# zymo <- readRDS("/Users/stephanieclouser/OneDrive - The Pennsylvania State University/Shared-Projects/1.Tutorials/ZymoAMR/zymo-amr.rds")

# get positive controls after removing contaminants
# pspos <- nocontam %>% 
#  ps_filter(str_detect(OperatorGender, "Pos")) 

# ntaxa(pspos)

# ntaxa(zymo)

# get_taxa_unique(pspos)
# get_taxa_unique(zymo)

# we don't want to remove the metals and biocides since we can't compare them to the zymo so we filter out only the drugs and multi-compounds
# pstruepos <- subset_taxa(pspos, Broadclass == "Drugs")
# get_taxa_unique(pstruepos)

# get what we want to remove
# psposrm <- subset_taxa(pstruepos, !taxa_names(pstruepos) %in% taxa_names(zymo))

# we need to remove the negative control from the sample set
pssave <- pssave <- ps_filter(nocontam, SampleBinary == "Sample")

# giving us the object decontaminated with all samples
pscount <- subset_taxa(ps, taxa_names(ps) %in% taxa_names(pssave))

# remove positive controls
ps.only.sample <- prune_samples(sample_data(pscount)$OperatorGender != "PositiveControl", pscount)
ps.only.sample2 <- prune_samples(sample_data(ps.only.sample)$OperatorGender != "NegativeControl", ps.only.sample)


# remove contaminant sequences
# fin_nocontam <- subset_taxa(ps.only.sample2, !taxa_names(pscount) %in% taxa_names(psposrm))

# removing SNP confirmation genes ----
psfilt <- ps.only.sample2 %>% tax_select(tax_list = "SNP", strict_matches = FALSE, deselect = TRUE)

Keep <- c("Male", "Female")

psfilt <- subset_samples(
  psfilt, 
  OperatorGender %in% Keep
)

# save the decontaminated phyloseq object for downstrem analysis
saveRDS(psfilt, file = "data/decontam-ps-fixed.rds")

# save work
save.image("data/decontam.RData")
