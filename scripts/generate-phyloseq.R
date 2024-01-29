#Integration of AMR++ output into phyloseq object for downstream analysis 
#Stephanie Bierly - 11/15/23 

#laod packges
library(BiocManager)
library(stringr)
library(tidyverse)
# install.packages("tidylog")
library(tidylog)
# BiocManager::install("phyloseq")
library(phyloseq)
# BiocManager::install("vegan", force = TRUE)
library(vegan) # version 2.6-4
library(ggplot2)
library(readxl)

#read in necessary files: count matrix, gene info, metadata---- 
countsDF <- read.delim("data/countmatrix-cleanedall.txt", sep = "\t") 
met <- read_excel("data/kenya-metadata.xlsx")
genes <- read.delim("data/geneinfo-all.txt", sep = "\t") %>%
  select(-c(MEG_ID)) %>%
  unique()

counts <- as.matrix(countsDF)
otutab <- phyloseq::otu_table(counts, taxa_are_rows = TRUE)

tax <- as.matrix(genes)
rownames(tax) <- genes$Gene
taxtab <- phyloseq::tax_table(tax)
taxa_names(taxtab)

samp <- phyloseq::sample_data(met)
rownames(samp) <- met$`SampleID`

ps <- phyloseq::phyloseq(otutab, taxtab, samp)


#save the ps object 
saveRDS(ps, "data/rawps.rds")
