# Set working directory
setwd("/Users/natal/R/whales")

library(ape)

# Read the CSV file (metadata)
labels <- read_csv("whale_metadata_minke.csv")
str(labels)

# Read the Newick tree
tree <- read.tree("2outgroup.treefile")

# Ensure tree tip labels are file names, not full paths
tree$tip.label <- basename(tree$tip.label)

# Prune to exclude minke
humpback_tree <- drop.tip(tree, tip = "humpback_sim_filt.bam")

# Root the tree with humpback as outgroup
humpback_tree <- root(humpback_tree, outgroup = "minke_sim_filt.bam", resolve.root = TRUE)

# Time calibration
# View nodes to identify the calibration point
N <- Ntip(humpback_tree)
root_node <- N + 1
print(paste("Root node:", root_node))  # Ensure node number for calibration is correct

# Create calibration data
calibration <- makeChronosCalib(humpback_tree, node = 76, age.min = 10.46, age.max = 10.50)

# Time-calibrate the tree
calibrated_tree <- chronos(humpback_tree, lambda = 1, calibration = calibration)
calibrated_tree

# Load required libraries
library("ggplot2")
library("ggtree")
library("treeio")
library("ape")
library("phytools")
library("tidyverse")

# Ensure the metadata matches the tree's tip labels
metadata <- labels %>%
  select(File, POP) %>%
  mutate(File = basename(File))  # Ensure File column matches tree tip labels

# Create a data frame for tree tips with metadata
tip_data <- data.frame(
  label = calibrated_tree$tip.label  # Extract tip labels from the tree
) %>%
  left_join(metadata, by = c("label" = "File"))  # Join metadata based on tip labels

# Calculate the maximum depth of the tree for proper scaling
max_time <- max(nodeHeights(calibrated_tree))

# Add metadata to the tree object and plot
p <- ggtree(calibrated_tree, layout = "rectangular") %<+% tip_data +  # Add metadata to the tree
  geom_tippoint(aes(color = POP), size = 3) +  # Replace tip labels with colored circles
  scale_color_manual(values = c("humpback" = "gray", "east" = "#1E88E5", "west" = "#e51e8a")) +  # Custom colors
  theme_tree2() +  # Add time scale
  labs(color = "Population") +  # Legend label
  ggtitle("Time Calibrated - Minke Outgroup") +  # Add title
  theme(legend.position = "right") +  # Adjust legend position
  scale_x_continuous(  # Keep the tree layout facing right
    limits = c(0, max_time),  # Set limits for 0 on the right and max_time on the left
    breaks = seq(0, max_time, length.out = 5),  # Adjust breaks
    labels = function(x) round(max_time - x, 2)  # Flip the time values
  ) +
  coord_cartesian(xlim = c(0, max_time))  # Ensure the plot area respects flipped axis

# Display the plot
print(p)
