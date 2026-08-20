# ibs analyses: nj tree, mds, clustering, distance distributions (Fig S8, 5D)
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")

library(ape)
library(ggplot2)
library(dplyr)
library(ggtree)

ids <- readLines("ids.txt")
M   <- as.matrix(read.table("whale_ibs.ibsMat"))
rownames(M) <- colnames(M) <- ids
pop <- read.csv("whale_pop_2026.csv", colClasses = c(sample_ID = "character"))          # sample_ID, pop
POP <- pop$pop[match(ids, pop$sample_ID)]
POP_colors <- c("east" = "#3FA5D5", "west" = "#E7AF23")

# ---- NJ tree (Fig S8) ----
tr <- nj(as.dist(M))
plot(tr, tip.color = POP_colors[POP], cex = 0.5, main = "IBS neighbour-joining tree")

# ---- PCA / MDS of IBS distances (Fig S8) ----
# ---- PCA / MDS of IBS distances (Fig S8) ----
mds <- cmdscale(as.dist(M), k = 2, eig = TRUE)

# Percent variance explained
var_explained <- round(
  100 * mds$eig[1:2] / sum(mds$eig[mds$eig > 0]),
  2
)

mdf <- data.frame(
  MDS1 = mds$points[,1],
  MDS2 = mds$points[,2],
  POP = POP
)

print(
  ggplot(mdf, aes(MDS1, MDS2, color = POP)) +
    geom_point(size = 3, alpha = 0.85) +
    scale_color_manual(values = POP_colors) +
    labs(
      title = "IBS MDS",
      x = paste0("PC1 (", var_explained[1], "%)"),
      y = paste0("PC2 (", var_explained[2], "%)"),
      color = "Sample Site"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    )
)

# ---- hierarchical clustering (Fig S8) ----
hc <- hclust(as.dist(M), method = "average")

# Convert to dendrogram
hc_tree <- as.phylo(hc)

p_hc <- ggtree(hc_tree, layout = "rectangular") %<+%
  data.frame(
    label = hc_tree$tip.label,
    POP = POP[match(hc_tree$tip.label, ids)]
  ) +
  geom_tippoint(
    aes(color = POP),
    size = 3
  ) +
  scale_color_manual(values = POP_colors) +
  theme_tree2() +
  theme(
    legend.position = "none"
  )

print(p_hc)

# ---- within / between-population IBS distance distributions (Fig 5D) ----
pairs <- t(combn(seq_along(ids), 2))
cat_lab <- apply(pairs, 1, function(p) {
  a <- POP[p[1]]; b <- POP[p[2]]
  if (a == b && a == "east") "within ENP"
  else if (a == b && a == "west") "within WNP"
  else "between"
})

dd <- data.frame(dist = M[pairs], comparison = cat_lab)

print(
  ggplot(dd, aes(dist, fill = comparison)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c(
      "within ENP" = "#3FA5D5",
      "within WNP" = "#E7AF23",
      "between"    = "grey70"
    )) +
    labs(
      x = "Genetic distance (1 - IBS)",
      y = "Density",
      title = "IBS distance distributions",
      fill = "Comparison"
    ) +
    theme_minimal() + 
    theme(
      axis.title = element_text(size = 16)
    )
)

