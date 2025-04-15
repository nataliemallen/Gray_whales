setwd("/Users/natal/R/whales")

library("ggplot2")     # For creating plots
library("ggtree")      # For visualizing phylogenetic trees
library("treeio")      # For handling phylogenetic tree data
library("reshape2")    # For reshaping data frames
library("ggstance")    # For horizontal dodging of ggplot elements
library("tidyverse")   # For data manipulation and visualization
library("ape")         # For phylogenetic analysis functions
library("phytools")


# Read the CSV file
labels <- read_csv("whale_metadata_2_outgroup.csv")
str(labels)

# Read the Newick tree
tree <- read.tree("2outgroup.treefile")

# Ensure tree labels are file names, not paths
tree$tip.label <- basename(tree$tip.label)
print(tree$tip.label)

# prune to exclude minke 
humpback_tree <- pruned.trees<-drop.tip(tree,tip="minke_sim_filt.bam") #exclude listed species
#write.tree(pruned.tree, file="pruned_phylogeny.tre") #save pruned tree to file

# make humpback outgroup
humpback_tree <- root(humpback_tree, outgroup = "humpback_sim_filt.bam", resolve.root = TRUE)
plot(humpback_tree)

# view nodes
str(humpback_tree)
N <- Ntip(humpback_tree)
root_node <- N + 1
root_node
# Plot the tree
plot(humpback_tree, show.tip.label = TRUE)
# Add node labels
nodelabels()

# time calibration
calibration <- makeChronosCalib(humpback_tree, node = 76, age.min = 7.47, age.max = 7.51)
calibrated_tree <- chronos(humpback_tree, lambda = 1, calibration = calibration)

# Plot the calibrated tree
plot(calibrated_tree, show.tip.label = TRUE, cex = 0.7)  # Adjust cex to resize labels
axisPhylo()  # Adds a time scale at the bottom






















str(calibrated_tree)
calibrated_tree$node.label <- paste0("Node_", seq_len(calibrated_tree$Nnode))
tree_tibble <- as_tibble(calibrated_tree)
colnames(labels)
colnames(tree_tibble)
tree1 <- full_join(tree_tibble, labels, by = c("label" = "File"))
# Convert the calibrated tree to a standard phylo object
phylo_tree <- as.phylo(calibrated_tree)

# Convert to treedata for ggtree
tree2 <- as.treedata(phylo_tree)

# Plot the tree
p1 <- ggtree(tree2, layout = "rectangular") +
  geom_tiplab(size = 3) +  # Add tip labels
  geom_treescale(x = 0, y = 0, fontsize = 3) +  # Add a scale bar
  geom_tippoint(aes(color = POP), size = 4, alpha = 1) +  # Color-coded tips
  scale_color_manual(values = c("east" = "#1E88E5", "west" = "#e51e8a", "humpback" = "gray")) +
  theme_tree2()

print(p1)



# Merge the calibrated tree with metadata
tree1 <- full_join(as_tibble(calibrated_tree), labels, by = c('label' = 'File'))

# Convert the merged tibble back into a phylogenetic tree object suitable for ggtree plotting
tree2 <- as.treedata(calibrated_tree)

# Create a time-scaled phylogenetic tree plot with color-coded tips
p1 <- ggtree(tree2, layout = "rectangular") +  # Use 'rectangular' layout for time scale compatibility
  geom_tiplab(size = 3) +                     # Add tip labels to the tree
  geom_treescale(x = 0, y = 0, fontsize = 3) +  # Add a scale bar
  geom_tippoint(aes(color = POP), size = 4, alpha = 1) +  # Add colored tip points based on 'POP'
  scale_color_manual(values = c("east" = "#1E88E5",       # Manually define colors
                                "west" = "#e51e8a",
                                "humpback" = "gray")) +
  theme_tree2()  # Use a theme with a time scale at the bottom

# Display the tree
print(p1)







tree <- root(tree, outgroup = "minke_sim_filt.bam", resolve.root = TRUE)
print(tree$edge.length)
##tree$edge.length[tree$edge.length == 0.0000000000] <- 1e-2

plot(tree)
print(tree)

# difference of 2.51
calibrated_tree <- chronos(tree, lambda = 1, calibration = makeChronosCalib(tree, node = 77, age.min = 10.46, age.max = 10.50))

plot(calibrated_tree)

# Plot the calibrated tree
plot(calibrated_tree, show.tip.label = TRUE, cex = 0.7)  # Adjust cex to resize labels
axisPhylo()  # Adds a time scale at the bottom







