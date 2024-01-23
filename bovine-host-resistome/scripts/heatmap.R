# Heatmap for ASM Abstract
# SAB 1/16/24

# load packages
library(phyloseq)
# install.packages("microbiome")
library(microbiome)
library(tibble)
library(dplyr)
# install.packages("OTUtable")
library(OTUtable)
library(viridis)
library(stringr)
library(readxl) # if your metadata is an excel sheet and not csv
library(pheatmap)
library(ggplot2)


# load phyloseq object
ps <- readRDS("bovine-host-resistome/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

# aggregate the taxa at the desired level (Broadclass, Class, Mechanism, or Gene)
psmech <- aggregate_taxa(psrel, level = "Class")

# export otu and tax tables to use for heatmap
otu <- as.data.frame(otu_table(psmech))
tax <- as.data.frame(tax_table(psmech))

# convert the rownames to a column to merge the tables
otu2 <- rownames_to_column(otu, "Class")

# merge the tax table and otu tables by common variable (desired level)
merge <- merge(otu2, tax, "Class")

# order the merged table so that the classes are sorted by Broadclass (for aesthetics of the plot so they can be grouped with separation between broadclasses)
m2 <- merge[order(merge$Broadclass),]

# remove unwanted columns, only select the samples and the class
df2 <- dplyr::select(m2, !c("Broadclass", "unique"))

# remove the rownames so that the class can be assigned to the rownames in the desired order
df2 <- rownames_to_column(df2, "random")

# remove the unwanted columns, only select the samples and the class
df2 <- dplyr::select(df2, !c("random"))

# assign classes to rownames
merge2 <- column_to_rownames(df2, "Class")

# calculate z scores of the OTU table we edited
z <- OTUtable::zscore(merge2)

# convert z scores to a matrix
z3 <- as.matrix(z)

# use the tax table to define the classes to broadclass
# select only the columns we want
merge3 <- dplyr::select(merge, c("Class", "Broadclass"))
# convert to data frame
merge3 <- as.data.frame(merge3)
# assign the class to rownames
merge3 <- column_to_rownames(merge3, "Class")
# rename the column "Broadclass" for aesthetic purposes on the heatmap
names(merge3)[names(merge3) == "Broadclass"] <- "Resistance Type"

#fix rownames on the matrix for aestheic purposes on the heatmap 
rownames(z3) <- str_replace_all(rownames(z3), "_", " ")#sub underscores with spaces
rownames(z3) <- str_replace(rownames(z3), "resistance", "")#drop resistance because that's a given
rownames(z3) <- str_replace_all(rownames(z3), "and", "+") #simplfy this to make shorter

# edit rownames in the broadclass info table to match z3 matrix
rownames(merge3) <- str_replace_all(rownames(merge3), "_", " ")#sub underscores with spaces
rownames(merge3) <- str_replace(rownames(merge3), "resistance", "")#drop resistance because that's a given
rownames(merge3) <- str_replace_all(rownames(merge3), "and", "+") #simplfy this to make shorter

# import metadata to assign variables at the top of the heatmap
# read in metadata
met <- read_excel("kenya-metadata.xlsx")
# convert the sample IDs to rownames
met <- column_to_rownames(met, "SampleID")
# select only the columns we want
met2 <- dplyr::select(met, c("OperatorGender", "HighestLevelEducation"))
# rename the columns for aesthetic purposes on the heatmap
names(met2)[names(met2) == "OperatorGender"] <- "Gender"
names(met2)[names(met2) == "HighestLevelEducation"] <- "Education"


# define colors for varables we want to see on the heatmap
colors <- list(
  "Gender" = c(Male = "#450c54", Female = "#24868EFF"),
  "Education" = c("Diploma or University degree" = "#eae69e", "Primary Level" = "#bfdb81", "Secondary Level" = "#83a561", "Tertiary education" = "#48723e"),
  "Resistance Type" = c("Biocides" = "#552F7A",
                        "Drugs" = "#7C5F98",
                        "Metals" = "#B09FC1",
                        "Multi-compound" = "#CABED6"))

# generate the heatmap
heat <- pheatmap(z3, 
         annotation_col = subset(x=met2, select = c("Gender", "Education")),
         annotation_row = subset(x=merge3, select = "Resistance Type"),
         # cluster_cols = hclust(dist(t(o6), method = "euclidean")), 
         cluster_rows = F,
         annotation_colors = colors,
         color = viridis(15),
         gaps_row = c(9,26,42),
         legend = T, 
        # fontsize = 15,
         annotation_legend = T,
         show_rownames = TRUE,
         border_color = NA)

# save the heatmap for export purposes
ggsave(heat, filename = "asm-heatmap-kenya.jpg", dpi = 600, width = 12, height = 12) # can save as .jpg or .pdf



