# manuscript phylogeny figure: rooted tree with bootstrap support and population tips
rm(list = ls())
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")

library(ggplot2)
library(ggtree)
library(treeio)
library(tidyverse)
library(ape)
library(phytools)

# ---- Read metadata and tree ----
labels <- read_csv("tree_tip_meta.csv", show_col_types = FALSE)

whale_tree <- read.tree("tree_final.treefile")
whale_tree$tip.label <- basename(whale_tree$tip.label)
whale_tree <- root(whale_tree, outgroup = "finwhale", resolve.root = TRUE)
whale_tree$edge.length[whale_tree$edge.length == 0] <- 1e-2

# ---- Defensive checks ----
message("Ntip: ", Ntip(whale_tree), "  Nnode: ", Nnode(whale_tree))
message("length(whale_tree$node.label): ", length(whale_tree$node.label))
message("First node labels (head):"); print(head(whale_tree$node.label, 10))

# Check metadata column that matches tip labels
if (!("File" %in% colnames(labels))) {
  stop("I expected the metadata tip-ID column 'File' to exist. If your tip IDs are in a different column, edit the script to use that column name for merging.")
}
# optional: ensure tip names match metadata
message("First few tip labels in tree:"); print(head(whale_tree$tip.label, 10))
message("First few 'File' values in metadata:"); print(head(labels$File, 10))

# ---- Fortify tree to table (coordinates + node ids + labels) ----
tree_df <- as_tibble(ggtree::fortify(whale_tree))

# fortify() produces a data.frame/tibble with columns including:
# node (numeric node id), label (tip labels and sometimes node labels), isTip, x, y, etc.

# ---- Add a node.label column matching internal nodes to whale_tree$node.label ----
nTip <- Ntip(whale_tree)
nNode <- Nnode(whale_tree)

# Initialize node.label NA
tree_df$node.label <- NA_character_

# In phylo, internal nodes are numbered from (nTip + 1) to (nTip + nNode).
# whale_tree$node.label is ordered for these internal nodes.
if (length(whale_tree$node.label) == nNode) {
  # For each row in tree_df that is internal, assign corresponding node.label
  internal_rows <- which(tree_df$node > nTip)
  # compute index into whale_tree$node.label:
  # for a node ID "node", index = node - nTip
  tree_df$node.label[internal_rows] <- whale_tree$node.label[tree_df$node[internal_rows] - nTip]
} else {
  warning("Number of node labels in whale_tree does not match Nnode(whale_tree). Node labels will be left NA.")
}

# ---- Merge metadata into tree_df by tip label ----
# We will join labels using: tree_df$label (tip names) == labels$File
# This will add POP, Path, etc to tip rows; internal nodes will get NA for these cols.
tree_df <- tree_df %>%
  left_join(labels, by = c("label" = "File"))

# Tip labels in visual top-to-bottom order
tip_order <- tree_df %>%
  filter(isTip) %>%
  arrange(desc(y)) %>%
  select(label, POP)

print(tip_order)

# Save if desired
write.csv(tip_order, "tree_tip_order_top_to_bottom.csv", row.names = FALSE)

# ---- Quick sanity checks (print a few rows) ----
message("Columns in fortified tree_df:"); print(colnames(tree_df))
message("First few tip rows (isTip == TRUE):")
print(tree_df %>% filter(isTip) %>% select(node, label, x, y, POP) %>% head())

message("First few internal node rows (isTip == FALSE):")
print(tree_df %>% filter(!isTip) %>% select(node, label, node.label, x, y) %>% head())

# ---- Build base ggtree plot (we will use tree_df for text/points so no %<+% needed) ----
base_plot <- ggtree(whale_tree, layout = "rectangular") +
  geom_treescale(x = 0, y = 60, offset.label = 50)

# find max x across tip x positions for the dotted extension
max_x <- max(tree_df$x[tree_df$isTip])

# ---- Final plot: add bootstrap labels, dotted tip extensions, and aligned tip circles ----
final_plot <- base_plot +
  # bootstrap labels: draw from tree_df internal-node rows
  geom_text(
    data = tree_df %>% 
      filter(!isTip & !is.na(node.label) & node.label != "Root"),
    aes(x = x, y = y, label = node.label),
    hjust = -0.2,
    size = 3
  ) +
  # dotted lines extending tips to max_x
  geom_segment(
    data = tree_df %>% filter(isTip),
    aes(x = x, xend = max_x, y = y, yend = y),
    linetype = "dotted", color = "black"
  ) +
  # single colored circle at end of dotted line per tip (POP from metadata)
  geom_point(
    data = tree_df %>% filter(isTip),
    aes(x = max_x, y = y, color = POP),
    size = 3, alpha = 1
  ) +

  scale_color_manual(
    values = c(
      "east" = "#3FA5D5",
      "west" = "#E7AF23",
      "fin" = "#3b3b3b"
    ),
    labels = c(
      "east" = "ENP",
      "west" = "WNP",
      "fin" = "Outgroup"
    )
  )+
  
  xlim(NA, max_x + max(tree_df$x) * 0.1)  # give some horizontal room for labels

# Print the plot
print(final_plot)

# --- Save high-resolution output files ---

# PDF output (PDF handles vector scaling automatically)
ggsave(
  filename = "whale_tree.pdf",
  plot = final_plot,
  width = 10, height = 12, units = "in"
)

# SVG output (dpi argument not relevant for pure vectors, but okay to include)
ggsave(
  filename = "whale_tree.svg",
  plot = final_plot,
  width = 10, height = 12, units = "in",
  dpi = 300
)

# get coordinates for admixture plots 
# Plot tree once to generate coordinates
plot(tr)

# Extract plotting coordinates
pp <- get("last_plot.phylo", envir = .PlotPhyloEnv)

# Get tip labels in top-to-bottom order
tip_order <- tr$tip.label[order(pp$yy[1:Ntip(tr)], decreasing = TRUE)]

print(tip_order)