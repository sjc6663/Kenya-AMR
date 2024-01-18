# load phyloseq object
ps <- readRDS("bovine-host-resistome/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

library(microbiome)
# WORKING CODE FOR HEATMAP ----
psmech <- aggregate_taxa(psrel, level = "Class")

otu <- as.data.frame(otu_table(psmech))
tax <- as.data.frame(tax_table(psmech))

library(tibble)
otu2 <- rownames_to_column(otu, "Class")

merge <- merge(otu2, tax, "Class")

m2 <- merge[order(merge$Broadclass),]

df2 <- dplyr::select(m2, !c("Broadclass", "unique"))

df2 <- rownames_to_column(df2, "random")
df2 <- dplyr::select(df2, !c("random"))
merge2 <- column_to_rownames(df2, "Class")

z <- OTUtable::zscore(merge2)

z3 <- as.matrix(z)

colors <- list(
  "Gender" = c(Male = "#450c54", Female = "#24868EFF"),
  "Education" = c("Diploma or University degree" = "#eae69e", "Primary Level" = "#bfdb81", "Secondary Level" = "#83a561", "Tertiary education" = "#48723e"),
  "Resistance Type" = c("Biocides" = "#552F7A",
                        "Drugs" = "#7C5F98",
                        "Metals" = "#B09FC1",
                        "Multi-compound" = "#CABED6"))

merge3 <- dplyr::select(merge, c("Class", "Broadclass"))
merge3 <- as.data.frame(merge3)
merge3 <- column_to_rownames(merge3, "Class")
names(merge3)[names(merge3) == "Broadclass"] <- "Resistance Type"

#fix rownames
library(stringr)
rownames(z3) <- str_replace_all(rownames(z3), "_", " ")#sub underscores with spaces
rownames(z3) <- str_replace(rownames(z3), "resistance", "")#drop resistance because that's a given
rownames(z3) <- str_replace_all(rownames(z3), "and", "+") #simplfy this to make shorter

# edit rownames to match z3 matrix
rownames(merge3) <- str_replace_all(rownames(merge3), "_", " ")#sub underscores with spaces
rownames(merge3) <- str_replace(rownames(merge3), "resistance", "")#drop resistance because that's a given
rownames(merge3) <- str_replace_all(rownames(merge3), "and", "+") #simplfy this to make shorter

library(readxl)
met <- read_excel("kenya-metadata.xlsx")
met <- column_to_rownames(met, "SampleID")
met2 <- dplyr::select(met, c("OperatorGender", "HighestLevelEducation"))
names(met2)[names(met2) == "OperatorGender"] <- "Gender"
names(met2)[names(met2) == "HighestLevelEducation"] <- "Education"


library(pheatmap)
library(viridis)
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

ggsave(heat, filename = "asm-heatmap-kenya.jpg", dpi = 600, width = 12, height = 12)



