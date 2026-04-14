# ==============================================================================
# MW Polyomavirus Phylogenetic Visualization
# Script: Phylogenetic Tree reconstruction and annotation
# Author: Thuany Giovana Daniel & Gabriel de Campos
# Requirements: ggtree, ape, phangorn, ggplot2, phytools, ggnewscale
# ==============================================================================

# Load Required Libraries
library(ggtree)      # Tree visualization
library(ggtext)      # Advanced text rendering
library(ape)         # Basic phylogenetics
library(ggplot2)     # Plotting framework
library(phangorn)    # Midpoint rooting
library(tidyr)       # Data manipulation
library(phytools)    # Phylogenetic tools
library(ggsci)       # Scientific color scales
library(ggnewscale)  # Multiple scales per plot

# Set Project Directory (Update to your project root)
# setwd("./phylogeny") 

# ------------------------------------------------------------------------------
# SECTION 1: Complete Genome Tree
# ------------------------------------------------------------------------------

# 1.1 Load and Root Tree
tree_complete <- read.tree("MWPyV_complete_consensus.fasta.treefile")
rooted_tree <- as.phylo(midpoint(tree_complete))

# 1.2 Load Metadata
metadata <- read.table("metadata_complete.tsv", sep="\t", header = TRUE)
metadata <- metadata[-1] # Remove index column if necessary

# 1.3 Core Plotting
# Attach metadata to tree structure
mw_plot <- ggtree(rooted_tree, color = "black") %<+% metadata 

# 1.4 Highlight Genotype Clades (Nodes based on IQ-TREE results)
mw_plot <- mw_plot + 
  geom_highlight(node = 40, fill = "#D4FCD6", alpha = 0.6) + # Genotype II
  geom_highlight(node = 38, fill = "#D4E2FC", alpha = 0.6) + # Genotype I
  geom_highlight(node = 59, fill = "#FCD9D4", alpha = 0.6)   # Genotype III

# 1.5 Final Aesthetic Layers
final_tree_complete <- mw_plot +
  # Sampling Location Color Palette
  scale_fill_manual(values = c("#CAB1BD", "#F0EC57", "#387780", "#221E36", 
                               "#E83151", "#D150B4", "#D2CCA1", "#BBCEA8"), 
                    name = "Sampling Location") +
  
  # Tip Labels (Study ID)
  geom_tiplab(aes(label = Study), size = 3.5, fontface = "bold.italic", 
              hjust = -0.15, offset = 0.0005) + 
  
  # Node Support (Bootstrap > 70%)
  geom_nodelab(
    aes(label = label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
    hjust = -0.15, size = 4, fontface = "bold"
  ) +
  
  # Map Country to Tip Points
  geom_tippoint(aes(fill = Country), size = 3, color = 'black', shape = 21, stroke = 0.5) +
  
  new_scale_fill() +
  geom_point(aes(x = 0, y = 0, fill = "Genotype I"), color = "transparent", size = 0) +
  geom_point(aes(x = 0, y = 0, fill = "Genotype II"), color = "transparent", size = 0) +
  geom_point(aes(x = 0, y = 0, fill = "Genotype III"), color = "transparent", size = 0) +
  
  scale_fill_manual(
    values = c("Genotype I" = "#D4E2FC", "Genotype II" = "#D4FCD6", "Genotype III" = "#FCD9D4"),
    name = "Genotypes",
    guide = guide_legend(order = 2, override.aes = list(color = "black", size = 5, shape = 22))
  ) +
  
  theme_tree2() +
  theme(legend.position = "right", plot.margin = margin(10, 120, 10, 20)) +
  expand_limits(x = 0.075)

# Save PDF (High Resolution)
ggsave("outputs/MWPyV_complete_genome.pdf", plot = final_tree_complete, width = 18, height = 9, dpi = 600)

# ------------------------------------------------------------------------------
# SECTION 2: Human Polyomavirus Comparison (MWPyV, JCPyV, BKPyV)
# ------------------------------------------------------------------------------

# 2.1 Load and Root
polyoma <- read.tree("HPyV_aln.fasta.treefile")
rooted_polyoma <- as.phylo(midpoint(polyoma)) 
metadados_hpyv <- read.table("metadata.tsv", sep="\t", header = TRUE)

# 2.2 Plotting Comparison Tree
p4 <- ggtree(rooted_polyoma, color = "lightgrey") %<+% metadados_hpyv +
  geom_tippoint(aes(fill = Virus), size = 3, color = 'black', shape = 21, stroke = 0.5) +
  scale_fill_manual(values = c("#1A659E", "#EFEFD0", "#FF6B35"), name = "Human Polyomavirus") +
  
  # Highlight target lineages (nodes 50 and 63)
  geom_tree(aes(color = (node %in% c(50, 63)))) +
  scale_color_manual(values = c("TRUE" = "#03191E", "FALSE" = "lightgrey")) +
  guides(color = "none") +
  
  geom_text(aes(label = Study), hjust = -.3, size = 3.5, fontface = "bold.italic") +
  geom_text2(aes(subset = node %in% c(50, 63), label = label), 
             hjust = -0.4, size = 4, fontface = "bold", color = "black") +
  
  theme_tree2() +
  theme(legend.position = "right", legend.text = element_text(face = "italic")) +
  expand_limits(x = 1.05)

# Save both PDF and JPEG
ggsave("outputs/HPyV_Tree.pdf", plot = p4, width = 18, height = 9, dpi = 600)

# ------------------------------------------------------------------------------
# SECTION 3: VP1 Gene Tree
# ------------------------------------------------------------------------------
# Note: Following the same logic as the Complete Genome but adjusted for VP1 nodes

tree_vp1 <- read.tree("VP1_MWPyV.fasta.treefile")
rooted_vp1 <- as.phylo(midpoint(tree_vp1))
metadata_vp1 <- read.table("metadata_VP1.tsv", sep="\t", header = TRUE)
metadata_vp1 <- metadata_vp1[-1]

final_tree_vp1 <- ggtree(rooted_vp1, color = "black") %<+% metadata_vp1 +
  geom_highlight(node = 61, fill = "#D4FCD6", alpha = 0.6) + # Genotype II
  geom_highlight(node = 37, fill = "#D4E2FC", alpha = 0.6) + # Genotype I
  geom_highlight(node = 69, fill = "#FCD9D4", alpha = 0.6) + # Genotype III
  scale_fill_manual(values = c("#CAB1BD", "#F0EC57", "#387780", "#221E36", 
                               "#E83151", "#D150B4", "#D2CCA1", "#BBCEA8"), 
                    name = "Sampling Location") +
  geom_tiplab(aes(label = Study), size = 3.5, fontface = "bold.italic", hjust = -0.45) +
  geom_nodelab(aes(label = label, subset = !is.na(as.numeric(label)) & as.numeric(label) > 70),
               hjust = -0.2, size = 4, fontface = "bold") +
  geom_tippoint(aes(fill = Country), size = 3, color = 'black', shape = 21, stroke = 0.5) +
  new_scale_fill() +
  scale_fill_manual(values = c("Genotype I" = "#D4E2FC", "Genotype II" = "#D4FCD6", "Genotype III" = "#FCD9D4"),
                    name = "Genotypes") +
  theme_tree2() +
  expand_limits(x = 0.065)

ggsave("outputs/MWPyV_VP1_Gene.pdf", plot = final_tree_vp1, width = 18, height = 9, dpi = 600)
