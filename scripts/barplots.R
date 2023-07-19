# Barplots for Descriptive Stats
# SAB 07/17/2023

# color scheme: Viridis mako / microshades micro_cvd_blue
color_palette <- c("#81007F", "#D75CE0", "#FF5EAA", "#FF8C76", "#FFC55A", "#F9F871")
short_color <- c("#D75CE0", "#FFC55A")

# load packages
library(phyloseq)
library(microViz)
library(ggplot2)
library(viridis)
library(BiocManager)
library(stringr)
library(microshades)
library(microbiome)

# read in phyloseq object
ps <- readRDS("data/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")

meltdf <- psmelt(psrel)

#filter out rows with 0 relabun 
out2 <- meltdf[meltdf$Abundance != 0.0000000000, ]

#plot 

# Male Female - Cows and Calves ----
ggplot(out2, aes(x=Broadclass, fill = HerdSize)) + 
  geom_bar()+
  theme_bw()+
  scale_fill_manual(values = color_palette) +
  facet_wrap(~OperatorGender, scales = "free", nrow = 2, ncol = 1)+
  labs(x = "",
       y = "Number of AMRg",
       fill = "Operator Gender",
       title = "C") +
  theme(axis.text= element_text(size = 12,face = "bold"),
        axis.title= element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 15, face = "bold", hjust = 0),
        legend.position = "right",
        strip.background = element_rect(fill = color_palette[5]),
        strip.text.x = element_text(size = 12, color = "black",face = "bold")
  ) +
  #guides(fill="none")+
  scale_y_continuous(expand = c(0,0.1), limits = c(0,5000))+
  coord_flip()

ggsave(filename = "plots/full-run/MF-calves-cows.pdf", dpi = 600)

# aggregate taxa ----
psbclass <- aggregate_taxa(psrel, level = "Broadclass")

# find and substitute
taxa_names(psbclass) <- gsub(taxa_names(psbclass), pattern = "_", replacement = " ") 

psbclass %>% plot_composition(group_by = "HighestLevelEducation", x.label = "OperatorGender") +
  # scale_y_continuous(labels = percent) +
  # theme(legend.position = "none") +
  scale_fill_manual(values = color_palette) +
  ggtitle("D")

ggsave(filename = "plots/full-run/relabund-averaged-bclass-gender.pdf", dpi = 600, width = 12, height = 14)
