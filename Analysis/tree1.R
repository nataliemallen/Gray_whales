# Load required libraries
library(ape)
library(phytools)
library(readr)
library(dplyr)

setwd("/Users/natal/R/whales")

# Read the CSV file
data <- read_csv("whale_metadata_5-24.csv")

# Read the Newick tree
tree <- read.tree("tree_3.nwk")

# Extract the file paths and pop information
file_paths <- data$Path
pop_info <- data %>% select(Path, POP)

# Create a named vector for the populations
pop_vector <- setNames(pop_info$POP, pop_info$Path)

# Replace file paths with population names in the tree's tip labels
tree$tip.label <- pop_vector[tree$tip.label]

# Manually set colors for each population using hex codes
pop_colors <- c("east" = "#1E88E5", "west" = "#D81B60")

# Assign colors to the tips based on population
tip_colors <- pop_colors[tree$tip.label]

# Plot the tree without tip labels
plot(tree, show.tip.label = FALSE)

# Add colored circles to the tips
tiplabels(pch = 21, bg = tip_colors, col = "black", cex = 1.5)

# Add legend
legend("topright", legend = names(pop_colors), fill = pop_colors, title = "POP")



###new tree with "east" and "west" tip labels 
# Create a named vector for the populations
pop_vector <- setNames(data$POP, data$Path)

# Replace the tip labels in the tree with the corresponding POP names
tree$tip.label <- pop_vector[tree$tip.label]

# Write the new tree to a Newick file
write.tree(tree, file="new_newicktree.nwk")


# Load required libraries
library(ape)
library(phytools)
library(readr)
library(dplyr)

setwd("/Users/natal/R/whales")

# Read the CSV file
data <- read_csv("whale_metadata_5-24.csv")

# Read the Newick tree
tree <- read.tree("tree_3.nwk")

# Extract the file paths and pop information
file_paths <- data$Path
pop_info <- data %>% select(Path, POP)

# Create a named vector for the populations
pop_vector <- setNames(pop_info$POP, pop_info$Path)

# Replace file paths with population names in the tree's tip labels
tree$tip.label <- pop_vector[tree$tip.label]

# Mid-root the tree
tree <- midpoint.root(tree)

# Manually set colors for each population using hex codes
pop_colors <- c("east" = "#1E88E5", "west" = "#D81B60")

# Assign colors to the tips based on population
tip_colors <- pop_colors[tree$tip.label]

# Plot the tree without tip labels
plot(tree, show.tip.label = FALSE)

# Add colored circles to the tips
tiplabels(pch = 21, bg = tip_colors, col = "black", cex = 1.5)

# Add legend
legend("topleft", legend = names(pop_colors), fill = pop_colors, title = "POP")
legend(x = 0, y = 0.95, legend = names(pop_colors), fill = pop_colors, title = "POP")

