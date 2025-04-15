# Set working directory
setwd("/Users/natal/R/whales")

library(ape)

# Read the CSV file (metadata)
labels <- read_csv("whale_metadata_2_outgroup.csv")
str(labels)

# Read the Newick tree
tree <- read.tree("2outgroup.treefile")

# Ensure tree tip labels are file names, not full paths
tree$tip.label <- basename(tree$tip.label)

# Prune to exclude minke
humpback_tree <- drop.tip(tree, tip = "minke_sim_filt.bam")

# Root the tree with humpback as outgroup
humpback_tree <- root(humpback_tree, outgroup = "humpback_sim_filt.bam", resolve.root = TRUE)

# Time calibration
# View nodes to identify the calibration point
N <- Ntip(humpback_tree)
root_node <- N + 1
print(paste("Root node:", root_node))  # Ensure node number for calibration is correct

# Create calibration data
calibration <- makeChronosCalib(humpback_tree, node = 76, age.min = 7.47, age.max = 7.51)

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
  ggtitle("Time Calibrated - Humpback Outgroup") +  # Add title
  theme(legend.position = "right") +  # Adjust legend position
  scale_x_continuous(  # Keep the tree layout facing right
    limits = c(0, max_time),  # Set limits for 0 on the right and max_time on the left
    breaks = seq(0, max_time, length.out = 5),  # Adjust breaks
    labels = function(x) round(max_time - x, 2)  # Flip the time values
  ) +
  coord_cartesian(xlim = c(0, max_time))  # Ensure the plot area respects flipped axis

# Display the plot
print(p)



###################################################################################################




# Ensure the tree and metadata are already loaded
# labels <- read_csv("whale_metadata_2_outgroup.csv")
# calibrated_tree <- chronos(humpback_tree, lambda = 1, calibration = calibration)

# Prepare metadata to match tree tip labels
metadata <- labels %>%
  select(File, POP) %>%
  mutate(File = basename(File))  # Ensure File column matches tree tip labels

# Create a data frame for tips
tip_data <- data.frame(
  label = calibrated_tree$tip.label  # Extract tip labels from the tree
) %>%
  left_join(metadata, by = c("label" = "File"))  # Join with metadata

# Plot the tree with color-coded tips
p <- ggtree(calibrated_tree) +
  geom_tiplab(aes(color = POP), size = 2, align = TRUE) +  # Add tip labels and color
  scale_color_manual(values = c("humpback" = "blue", "east" = "green", "west" = "red")) +  # Custom colors
  theme_tree2() +  # Add time scale
  labs(color = "Population") +  # Legend label
  ggtitle("Time-Calibrated Phylogenetic Tree") +  # Add title
  theme(legend.position = "right")  # Adjust legend position

# Display the plot
print(p)


# Ensure the labels file and tree are already loaded
# labels <- read_csv("whale_metadata_2_outgroup.csv")
# calibrated_tree <- chronos(humpback_tree, lambda = 1, calibration = calibration)

# Convert the tree to a tibble, focusing only on the tips
tree_tibble <- as_tibble(calibrated_tree)

# Filter only the tip labels (these correspond to the leaves of the tree)
tip_data <- tree_tibble %>%
  filter(isTip) %>%
  select(label) %>%
  left_join(metadata, by = c("label" = "File"))  # Join metadata based on tip labels

# Plot the tree with color-coded tips
p <- ggtree(calibrated_tree) %<+% tip_data +  # Add metadata to the tree plot
  geom_tiplab(aes(color = POP), size = 2, align = TRUE) +  # Add tip labels and color
  scale_color_manual(values = c("humpback" = "blue", "east" = "green", "west" = "red")) +  # Define colors
  theme_tree2() +  # Add time scale
  labs(color = "Population") +  # Legend label
  ggtitle("Time-Calibrated Phylogenetic Tree") +  # Add title
  theme(legend.position = "right")  # Adjust legend position

# Display the plot
print(p)


# Create a data frame mapping tip labels to metadata
metadata <- labels %>%
  select(File, POP) %>%
  mutate(File = basename(File))  # Ensure the File column matches tree tip labels

# Merge metadata into the tree data
tree_data <- as_tibble(calibrated_tree) %>%
  left_join(metadata, by = c("label" = "File"))

# Plot the tree with a time scale and color tips by "POP"
ggtree(calibrated_tree, aes(color = POP)) %>%
  ggtree::geom_tiplab(size = 2, align = TRUE) +  # Add tip labels
  scale_color_manual(values = c("humpback" = "blue", "east" = "green", "west" = "red")) +  # Custom colors
  theme_tree2() +  # Time scale theme
  labs(color = "Population") +  # Legend label
  ggtitle("Time-Calibrated Phylogenetic Tree") +  # Add title
  theme(legend.position = "right")  # Adjust legend position
























tree_phylo <- as.phylo(calibrated_tree)

# Join metadata with tree tip labels
tip_metadata <- data.frame(label = calibrated_tree$tip.label) %>%
  left_join(labels, by = c("label" = "File"))

print(head(tip_metadata))

library(ggtree)

# Plot the tree
p <- ggtree(tree_phylo, layout = "rectangular") +
  # Add a time scale
  geom_treescale(x = 0, y = 0, fontsize = 3) +
  # Add colored tip points based on population (POP)
  geom_tippoint(aes(color = POP), data = tip_metadata, size = 4, alpha = 1) +
  # Define colors for populations
  scale_color_manual(values = c(
    "east" = "#1E88E5",   # Blue for east
    "west" = "#e51e8a",  # Magenta for west
    "humpback" = "gray"  # Gray for humpback
  )) +
  # Customize tree theme
  theme_tree2()

# Print the plot
print(p)











# Join metadata with tree tip labels
tip_metadata <- data.frame(label = calibrated_tree$tip.label) %>%
  left_join(labels, by = c("label" = "File"))
# Ensure metadata matches the tree tip labels
filtered_labels <- labels %>%
  filter(File %in% calibrated_tree$tip.label)  # Filter metadata to match tree tips

# Add a `POP` mapping to tips
tip_metadata <- data.frame(label = calibrated_tree$tip.label) %>%
  left_join(filtered_labels, by = c("label" = "File"))

# Plot the tree
p <- ggtree(calibrated_tree, layout = "rectangular") +
  # Add a time scale at the bottom
  geom_treescale(x = 0, y = 0, fontsize = 3) +
  # Add tip points with metadata-based colors
  geom_tippoint(aes(color = POP), data = tip_metadata, size = 4, alpha = 1) +
  # Define the colors for each population
  scale_color_manual(values = c("east" = "#1E88E5",   # Blue for east
                                "west" = "#e51e8a",  # Magenta for west
                                "humpback" = "gray")) +  # Gray for humpback
  # Add tree theme
  theme_tree2()

# Print the plot
print(p)


# Ensure the metadata matches the tree tip labels
labels <- labels %>%
  filter(File %in% calibrated_tree$tip.label)  # Only keep matching tips

# Add metadata to the tree for plotting
calibrated_tree$tip.label <- as.character(calibrated_tree$tip.label)  # Ensure tip labels are character

# Plot the tree directly
p <- ggtree(calibrated_tree, layout = "rectangular") +
  # Add a time scale at the bottom
  geom_treescale(x = 0, y = 0, fontsize = 3) +
  # Add tip points colored by population (POP)
  geom_tippoint(aes(subset = calibrated_tree$tip.label %in% labels$File, 
                    color = labels$POP[match(calibrated_tree$tip.label, labels$File)]), 
                size = 4, alpha = 1) +
  # Manually define the colors for each population
  scale_color_manual(values = c("east" = "#1E88E5",   # Blue for east
                                "west" = "#e51e8a",  # Magenta for west
                                "humpback" = "gray")) +  # Gray for humpback
  # Add tree theme
  theme_tree2()

# Print the plot
print(p)

# Check structure of labels
str(labels)

# Print tip labels of the calibrated tree
print(calibrated_tree$tip.label)

# Check if all tree tip labels exist in the metadata
missing_labels <- setdiff(calibrated_tree$tip.label, labels$File)
missing_labels
if (length(missing_labels) > 0) {
  print("Missing labels in metadata:")
  print(missing_labels)
}

print(calibrated_tree)
