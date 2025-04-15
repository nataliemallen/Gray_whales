# Load required libraries
library(ape)
library(phytools)
library(readr)
library(dplyr)

setwd("/Users/natal/R/whales")


# unrooted tree -----------------------------------------------------------

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

#####nuclear tree (mid-rooted)
# Set the working directory to the folder containing your R project files
setwd("/Users/natal/R/whales")

# Load necessary libraries for plotting, tree manipulation, and data wrangling
library("ggplot2")     # For creating plots
library("ggtree")      # For visualizing phylogenetic trees
library("treeio")      # For handling phylogenetic tree data
library("reshape2")    # For reshaping data frames
library("ggstance")    # For horizontal dodging of ggplot elements
library("tidyverse")   # For data manipulation and visualization
library("ape")         # For phylogenetic analysis functions
library("phytools")

# Read the CSV file
labels <- read_csv("whale_metadata_5-24.csv")
#View(labels)

# Read the Newick tree
turtle_tree <- read.tree("tree_3.nwk")

# Create a named vector where the names are the file paths and the values are the corresponding Fnumbers
path_to_fnumber <- setNames(labels$File, labels$Path)

# Replace the file paths in the tree with the corresponding Fnumbers
turtle_tree$tip.label <- path_to_fnumber[turtle_tree$tip.label]

# Midroot the tree at its midpoint
turtle_tree <- midpoint.root(turtle_tree)

# Merge the updated tree data (converted to tibble) with the labels dataframe using 'label' in the tree
# and 'Fnumber' in the labels as the key columns
tree1 <- full_join(as_tibble(turtle_tree), labels, by = c('label' = 'File'))

# Convert the merged tibble back into a phylogenetic tree object suitable for ggtree plotting
tree2 <- as.treedata(tree1)

# Create a circular phylogenetic tree plot using the tree2 object
p1 <- ggtree(tree2, layout = "circular", branch.length = "none") +
  # Add a scale bar to the tree at the specified coordinates
  geom_treescale(x=0, y=60, offset.label = 50) +
  # Add colored tip points to the tree, colored by the 'Site' variable
  geom_tippoint(aes(color = POP), size = 4, alpha = 1, position = "identity") +
  # Manually define the colors associated with each site name
  scale_color_manual(values = c("east" = "#1E88E5",  
                                "west" = "#e51e8a"))

#p1 <- ggtree(tree2, layout = "circular", branch.length = "none") +
#  # Add a scale bar to the tree at the specified coordinates
#  geom_treescale(x=0, y=60, offset.label = 50) +
#  # Add colored tip points to the tree, colored by the 'Site' variable
#  geom_tippoint(aes(color = River), size = 4, alpha = 1, position = "identity") +
#  # Manually define the colors associated with each site name
#  scale_color_manual(values = c("Neches" = "#1E88E5",  # Deep blue
#                                "Brazos" = "#56B4E9",       # Light blue
#                                "San Jacinto" = "#009E73",       # Green
#                                "Trinity" = "#8E44AD",        # Purple
#                                "Colorado" = "#d13b68") )

# Display the final tree plot
p1


#####nuclear tree with outgroup
# Set the working directory to the folder containing your R project files
setwd("/Users/natal/R/whales")

# Load necessary libraries for plotting, tree manipulation, and data wrangling
library("ggplot2")     # For creating plots
library("ggtree")      # For visualizing phylogenetic trees
library("treeio")      # For handling phylogenetic tree data
library("reshape2")    # For reshaping data frames
library("ggstance")    # For horizontal dodging of ggplot elements
library("tidyverse")   # For data manipulation and visualization
library("ape")         # For phylogenetic analysis functions
library("phytools")

# Read the CSV file
labels <- read_csv("whale_metadata_5-24_outgroup.csv")
#View(labels)

# Read the Newick tree
turtle_tree <- read.tree("outgroup.not.iupac.resolved.treefile")

turtle_tree$tip.label <- basename(turtle_tree$tip.label)
print(turtle_tree$tip.label)

# Create a named vector where the names are the file paths and the values are the corresponding Fnumbers
#path_to_fnumber <- setNames(labels$File, labels$Path)

# Replace the file paths in the tree with the corresponding Fnumbers
#turtle_tree$tip.label <- path_to_fnumber[turtle_tree$tip.label]
#print(turtle_tree$tip.label)

# Midroot the tree at its midpoint
#turtle_tree <- midpoint.root(turtle_tree)

# Merge the updated tree data (converted to tibble) with the labels dataframe using 'label' in the tree
# and 'Fnumber' in the labels as the key columns
tree1 <- full_join(as_tibble(turtle_tree), labels, by = c('label' = 'File'))

# Convert the merged tibble back into a phylogenetic tree object suitable for ggtree plotting
tree2 <- as.treedata(tree1)

# Create a circular phylogenetic tree plot using the tree2 object
p1 <- ggtree(tree2, layout = "circular", branch.length = "none") +
  # Add a scale bar to the tree at the specified coordinates
  geom_treescale(x=0, y=60, offset.label = 50) +
  # Add colored tip points to the tree, colored by the 'Site' variable
  geom_tippoint(aes(color = POP), size = 4, alpha = 1, position = "identity") +
  # Manually define the colors associated with each site name
  scale_color_manual(values = c("east" = "#1E88E5",  
                                "west" = "#e51e8a",
                                "fin" = "gray"))

# Display the final tree plot
p1


#### new
# Set the working directory to the folder containing your R project files
setwd("/Users/natal/R/whales")

# Load necessary libraries for plotting, tree manipulation, and data wrangling
library("ggplot2")     # For creating plots
library("ggtree")      # For visualizing phylogenetic trees
library("treeio")      # For handling phylogenetic tree data
library("reshape2")    # For reshaping data frames
library("ggstance")    # For horizontal dodging of ggplot elements
library("tidyverse")   # For data manipulation and visualization
library("ape")         # For phylogenetic analysis functions
library("phytools")    # For phylogenetic manipulation

# Read the CSV file
labels <- read_csv("whale_metadata_5-24_outgroup.csv")

# Read the Newick tree
turtle_tree <- read.tree("outgroup.not.iupac.resolved.treefile")

# Ensure tree labels are file names, not paths
turtle_tree$tip.label <- basename(turtle_tree$tip.label)
print(turtle_tree$tip.label)

turtle_tree <- root(turtle_tree, outgroup = "SRR23615109_filt.bam", resolve.root = TRUE)
print(turtle_tree$edge.length)
turtle_tree$edge.length[turtle_tree$edge.length == 0.0000000000] <- 1e-2

#tree
ggtree(turtle_tree, layout = "rectangular") %<+% labels +
  geom_treescale(x = 0, y = 60, offset.label = 50) +
  geom_tippoint(aes(color = POP), size = 3, alpha = 1) +
  scale_color_manual(values = c("east" = "#1E88E5",
                                "west" = "#e51e8a",
                                "fin" = "darkgray"))

# to get tree with dotted lines extending leaves
p <- ggtree(turtle_tree, layout = "rectangular") %<+% labels + 
  geom_treescale(x = 0, y = 60, offset.label = 50) + 
  geom_tippoint(aes(color = POP), size = 3, alpha = 1) + 
  scale_color_manual(values = c("east" = "#1E88E5", "west" = "#e51e8a", "fin" = "darkgray"))

# extract the tree data from the ggtree object
tree_data <- p$data

# find the maximum x position (tree length)
max_x <- max(tree_data$x[tree_data$isTip])

# add segments to extend tips
p + geom_segment(data = tree_data[tree_data$isTip, ], aes(
  x = x,
  xend = max_x,
  y = y,
  yend = y
), linetype = "dotted", color = "black")

#### tree with circles lined up 
# get the ggtree object
p <- ggtree(turtle_tree, layout = "rectangular") %<+% labels + 
  geom_treescale(x = 0, y = 60, offset.label = 50)

# extract the tree data from the ggtree object
tree_data <- p$data

# find the maximum x position (tree length)
max_x <- max(tree_data$x[tree_data$isTip])

# add segments to extend tips ### use this tree
p + 
  geom_segment(data = tree_data[tree_data$isTip, ], aes(
    x = x,
    xend = max_x,
    y = y,
    yend = y
  ), linetype = "dotted", color = "black") +
  # move the colored circles to the end of the dashed line
  geom_point(data = tree_data[tree_data$isTip, ], aes(
    x = max_x,
    y = y,
    color = POP
  ), size = 3, alpha = 1) +
  # ensure color scaling matches the labels dataset
  scale_color_manual(values = c("east" = "#1E88E5", "west" = "#e51e8a", "fin" = "darkgray"))

#####
# get the ggtree object
p <- ggtree(turtle_tree, layout = "rectangular") %<+% labels + 
  geom_treescale(x = 0, y = 60, offset.label = 50)

# extract the tree data from the ggtree object
tree_data <- p$data

# find the maximum x position (tree length)
max_x <- max(tree_data$x[tree_data$isTip])

# add segments to extend tips
p + 
  geom_segment(data = tree_data[tree_data$isTip, ], aes(
    x = x,
    xend = max_x,
    y = y,
    yend = y
  ), linetype = "dotted", color = "black") +
  # move the colored circles to the end of the dashed line
  geom_point(data = tree_data[tree_data$isTip, ], aes(
    x = max_x,
    y = y,
    color = POP
  ), size = 3, alpha = 1) +
  # ensure color scaling matches the labels dataset
  scale_color_manual(values = c("east" = "#1E88E5", "west" = "#e51e8a", "fin" = "darkgray")) +
  # add the tip labels to the right of the colored circles
  geom_text(data = tree_data[tree_data$isTip, ], aes(
    x = max_x - 0.05,  # increase the space for labels
    y = y,
    label = label
  ), hjust = 0, size = 3) +
  # move the legend to the upper left corner
  theme(legend.position = c(0.05, 0.95))

