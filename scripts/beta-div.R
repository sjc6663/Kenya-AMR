# beta diversity 
# SAB 07/17/2023

# load packages
library(phyloseq)
library(microViz)
library(vegan)
library(ggplot2)
library(patchwork)

set.seed(81299)

color_palette <- c("#81007F", "#D75CE0", "#FF5EAA", "#FF8C76", "#FFC55A", "#F9F871")
short_color <- c("#D75CE0", "#FFC55A")

# load phyloseq object
ps <- readRDS("data/decontam-ps.RDS")

# transform to relative abundance
psrel <- microbiome::transform(ps, "compositional")


# clr transform phyloseq objects
transps <- psrel %>% 
  tax_transform(trans = "clr") %>% 
  ps_get()


dist_mat <- phyloseq::distance(transps, method = "euclidean")

vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$OperatorGender) # p = 0.845, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$HerdSize) # p = 0.076, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$NonFamilyEmployees) # p = 0.83, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$LastABXUse) # p = 0.68, ns
vegan::adonis2(dist_mat ~ phyloseq::sample_data(transps)$HighestLevelEducation) # p = 0.58, ns

# PCA Plots ------------------------------------------------------------------------------------------------------

# Gender
A <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "OperatorGender") +
  scale_color_manual(values = short_color) +
  stat_ellipse(aes(group = OperatorGender, color = OperatorGender)) + 
  theme_classic() +
  ggtitle("A") 
  #labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")

# Non-Family Employees
B <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "NonFamilyEmployees") +
  scale_color_manual(values = short_color) +
  stat_ellipse(aes(group = NonFamilyEmployees, color = NonFamilyEmployees)) + 
  theme_classic() +
  ggtitle("B") 
  #labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")

# Herd Size
C <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "HerdSize") +
  scale_color_manual(values = color_palette) +
  stat_ellipse(aes(group = HerdSize, color = HerdSize)) + 
  theme_classic() +
  ggtitle("C") 
  #labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")

# Last ABX Use
D <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "LastABXUse") +
  scale_color_manual(values = color_palette) +
  stat_ellipse(aes(group = LastABXUse, color = LastABXUse)) + 
  theme_classic() +
  ggtitle("D")  
  #labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")

# Highest Level Edu
E <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "HighestLevelEducation") +
  scale_color_manual(values = color_palette) +
  stat_ellipse(aes(group = HighestLevelEducation, color = HighestLevelEducation)) + 
  theme_classic() +
  ggtitle("E") 
  #labs(caption = "R2 = 0.016, F(1,70) = 1.07, P = 0.30")

(A/B/C)|(D/E)

ggsave(filename = "plots/beta-div-all.pdf", dpi = 600, width = 15, height = 16)

# beta diversity ----

ait <- ps %>%
  # transform to relative abundance
  tax_transform("compositional", keep_counts = FALSE) %>%
  dist_calc("aitchison")


# test beta dispersion
ait %>% dist_bdisp(variables = "OperatorGender") %>% bdisp_get() # p=0.105

# test with PERMANOVA
mod1 <- ait %>%
  dist_permanova(
    seed = 81299,
    variables = "OperatorGender",
    n_perms = 9999
  )

mod1 # R2 = 0.03, F(1, 65) = 1.88, P = 0.0001***

C1 <- ait %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method, in this case we are using the already transformed object above with aitchison's distance, so we use PCoA
  # tax_transform(trans = "clr") %>%
  ord_calc(method = "PCoA") %>% 
  ord_plot(color = "OperatorGender", size = 6, axes = c(2,3)) +
  scale_color_manual(values = c("#D75CE0", "#FFC55A")) +
  stat_ellipse(aes(group = OperatorGender, color = OperatorGender)) + 
  theme_classic() +
  labs(color = "Gender") +
  ggtitle("C") + 
  #labs(caption = "R2 = 0.03, F(1, 65) = 1.88, P = 0.0001***") +
  theme(text = element_text(size = 20)) 
C1

C2 <- psrel %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  tax_transform(trans = "clr") %>%
  ord_calc(method = "PCA") %>% 
  ord_plot(color = "OperatorGender", size = 6, axes = c(2,3)) +
  scale_color_manual(values = c("#D75CE0", "#FFC55A")) +
  stat_ellipse(aes(group = OperatorGender, color = OperatorGender)) + 
  theme_classic() +
  labs(color = "Gender") +
  ggtitle("C") + 
 # labs(caption = "R2 = 0.03, F(1, 65) = 1.88, P = 0.0001***") +
  theme(text = element_text(size = 20)) 
C2

library(patchwork)
C1|C2
# this shows that these different codes produce essentially the same plot, which makes me feel a lot better about what I've been doing, lol.