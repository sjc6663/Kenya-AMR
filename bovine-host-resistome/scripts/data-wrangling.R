## Gene Info Table
## SAB last updated 07/11/2023

#load packages----
library(BiocManager)
library(dplyr)
library(stringr)

#read in and clean data----
b1 <- read.csv("amr++outputs/bovine-host/AMR_analytic_matrix.csv")

#do this for all batches, n being the batch n 
df1 <- cbind(b1, as.data.frame(str_split_fixed(rownames(b1), "\\|", 5))[5]) %>%
  unique() %>%
  aggregate(.~V5, FUN= sum)

# install.packages("plyr")
library(plyr)
merged <- df1
merged[is.na(merged)] <- 0 #replace NAs 
rownames(merged) <- merged$V5 #fix rownames
colnames(merged) <- trimws(colnames(merged), whitespace = "_.*") #fix colnames
merged_clean <- merged[2:ncol(merged)]

#write.table(merged_clean, file = "bovine-host-resistome/countmatrix-cleanedall.txt", sep = "\t") #save this matrix


#break up gene information for alllll batches - this will be useful for later 
#install.packages("rlist")
library("rlist")
genes <- as.data.frame(list.append(rownames(b1))) %>%
  unique()

colnames(genes)[1] <- "allmeg"
genes <- as.data.frame(str_split_fixed(genes$allmeg, "\\|", 5))
colnames(genes) <- c("MEG_ID", "Broadclass", "Class", "Mechanism", "Gene")

#write.table(genes, file = "bovine-host-resistome/geneinfo-all.txt", sep = "\t") 
