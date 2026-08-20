# pca plots from pcangsd covariance matrices (Fig 4B, S2-S4)
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")

library(ggplot2)
library(dplyr)

# Metadata MUST be in the same order as bams_list.txt (73 individuals, 038046 removed).
# Columns expected: File, POP, Order  (POP = "east"/"west")
meta <- read.csv("whale_metadata_2026.csv")

POP_colors <- c("east" = "#3FA5D5", "west" = "#E7AF23")

# Helper: read a PCAngsd .cov, run prcomp, attach metadata, return df + var explained
pca_from_cov <- function(cov_file, meta_df) {
  cov <- as.matrix(read.table(cov_file))
  pr  <- prcomp(cov, scale = TRUE)
  ve  <- (pr$sdev^2 / sum(pr$sdev^2)) * 100
  df  <- cbind(meta_df, as.data.frame(pr$x))
  list(df = df, ve = ve)
}

plot_pca <- function(res, title) {
  ggplot(res$df, aes(x = PC1, y = PC2, color = POP)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = POP_colors) +
    labs(x = paste0("PC1 (", round(res$ve[1], 2), "%)"),
         y = paste0("PC2 (", round(res$ve[2], 2), "%)"),
         color = "Sample Site", title = title) +
    theme_minimal() +
    theme(axis.title = element_text(size = 14), legend.title = element_text(size = 12))
}

# ---- Fig 4B: genome-wide ----
gw <- pca_from_cov("whale.cov", meta)
print(plot_pca(gw, "PCA - all individuals"))

# ---- Fig S2: LD-pruned ----
pr <- pca_from_cov("whale_pruned.cov", meta)
print(plot_pca(pr, "PCA - LD pruned"))

# ---- Fig S2: relatives excluded ----
# meta_norel must match norel_bams.txt order
remove     <- readLines("relatives_remove.txt")
meta_norel <- meta[!sub("\\.bam$", "", meta$File) %in% remove, ]
write.csv(meta_norel, "whale_metadata_2026_norel.csv", row.names = FALSE)

meta_norel <- read.csv("whale_metadata_2026_norel.csv")
nr <- pca_from_cov("whale_norel.cov", meta_norel)
print(plot_pca(nr, "PCA - relatives excluded"))

# ---- Fig S3: per sample site ----
meta_east <- subset(meta, POP == "east")
meta_west <- subset(meta, POP == "west")
print(
  plot_pca(pca_from_cov("east.cov", meta_east), "PCA - ENP only") +
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5))
)

print(
  plot_pca(pca_from_cov("west.cov", meta_west), "PCA - WNP only") +
    coord_cartesian(xlim = c(-5, 5), ylim = c(-5, 5))
)
# ---- Fig S4: per chromosome ----
outdir <- "chr_pca_plots"
dir.create(outdir, showWarnings = FALSE)

chrs <- readLines("chrs.txt")

for (chr in chrs) {
  f <- paste0("chr_", chr, ".cov")
  
  if (file.exists(f)) {
    
    p <- plot_pca(
      pca_from_cov(f, meta),
      paste0("PCA - ", chr)
    ) +
      labs(x = NULL, y = NULL) +          # Remove PC1 and PC2 labels
      theme(
        plot.title = element_text(size = 20, hjust = 0.5),
        legend.position = "none"
      )
    
    ggsave(
      filename = file.path(outdir, paste0("PCA_", chr, ".png")),
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
}
