# Kenya Dairy Resistome 
Analysis of whole shotgun metagenomics data from samples collected on dairy farms in Kenya. This project desires to look at the human and ethical dimensions of antimicrobial resistance (AMR) in Kenya. This repository includes scripts related to analyzing that data. Analysis was done using AMR++ for gene alignment and downstream analysis was completed with R in RStudio. Below is a breakdown of what each file is relevant for related to analysis.

This data is currently being analyzed and written for publication. 

## Data

### countmatrix-cleanedall.txt
This file contains the count information for all AMRg aligned via the AMRPlusPlus pipeline. 

### geneinfo-all.txt
This file contains the gene information for all AMRg aligned via the AMRPlusPlus pipeline. 

### kenya-metadata.xlsx
This file contains the anonymous metadata for all samples in the study. 

### decontam-ps.rds
This file is the post-decontam phyloseq object used for all analysis of data. 

## Scripts

### data-wrangling.R
This R script takes the output from AMR++ and separates the gene info and the count matrix info that are used to generate the phyloseq object for downstream analysis. 

### generate-phyloseq.R
This R script takes the count matrix, the gene info, and the metadata and combines it into a phyloseq object that is used for downstream analyses. 

### decontam.R
This R script utilizes the decontam package to remove putative contaminant AMRg found in the negative controls sequenced with samples. 

### gender.R
This R script runs statistical tests for alpha and beta diversity by gender and creates a relative abundance plot, a boxplot for alpha diversity, and a PCoA for beta diversity, all for gender of primary farm operator. 

### education.R
This R script runs statistical tests for alpha and beta diversity by education level and creates a relative abundance plot, a boxplot for alpha diversity, and a PCoA for beta diversity, all for highest level of education of primary farm operator. 

### herdsize-3.R
This R script runs statistical tests for alpha and beta diversity by herd size (greater than or less than/equal to 3) and creates a relative abundance plot, a boxplot for alpha diversity, and a PCoA for beta diversity, all for herd size. 

### heatmap.R
This R script generates a heatmap for AMRg at the Class level for all samples. It also adds identifier information to the samples for gender and highest level of education for primary farm operator. 

### presence-absence-binary.R
This R script converts the OTU table to a binary format to use for presence/absence of AMRg for analysis. 
