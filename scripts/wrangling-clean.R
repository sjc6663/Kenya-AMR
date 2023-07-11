## Wrangling and Cleanup AMR++ Output
## SAB 07/11/2023

# load packages
library(BiocManager)
library(dplyr)
library(stringr)

# read in and clean data
out <- read.csv("kenya_723_AMR_analytic_matrix.csv")

# drop extension from file name
colnames(out) <- str_split_i(colnames(out), "_", 1)


# break up gene information
genes <- as.data.frame(rownames(out))
genes <- as.data.frame(str_split_fixed(genes$`rownames(out)`, "\\|", 5))
# give the column names
colnames(genes) <- c("MEG_ID", "Broadclass", "Class", "Mechanism", "Gene")

out2 <- cbind(out, genes[5])

# collapse the genes so that they are only one row
out3 <- aggregate(.~Gene, data = out2, FUN = sum)
# make the table super pretty
out4 <- out3[2:ncol(out3)]
# change the rownames to be the genes
rownames(out4) <- out3$Gene
# save the file 
write.table(out4, file = "data/countmatrix-cleaned.txt", sep = "\t")
