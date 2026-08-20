# sliding-window pca: per-window pc1 along each chromosome (Fig S5)

setwd("/scratch/gautschi/allen715/2026_whales/sliding_pca")

library(ggplot2)
library(dplyr)

POP_colors <- c("east" = "#3FA5D5", "west" = "#E7AF23")

# Output directory
outdir <- "sliding_pca_plots"
dir.create(outdir, showWarnings = FALSE)

# meta must match whale.beagle.gz individual order (73 individuals)
meta <- read.csv("whale_metadata_2026.csv")
man  <- read.table(
  "window_manifest.txt",
  col.names = c("chr", "start", "nsnp", "cov")
)

# ------------------------------------------------------------------
# Extract PC1 (per individual) for every window
# ------------------------------------------------------------------

recs <- vector("list", nrow(man))

for (i in seq_len(nrow(man))) {

  f <- man$cov[i]

  if (!file.exists(f))
    next

  cov <- as.matrix(read.table(f))
  e <- eigen(cov)
  pc1 <- e$vectors[, 1]

  recs[[i]] <- data.frame(
    chr = man$chr[i],
    pos = man$start[i] + 50000,   # window midpoint
    ind = seq_along(pc1),
    POP = meta$POP,
    PC1 = pc1
  )
}

dat <- bind_rows(recs)

# ------------------------------------------------------------------
# One plot per chromosome with east/west facets
# ------------------------------------------------------------------

for (cc in unique(dat$chr)) {

  d <- subset(dat, chr == cc)

  p <- ggplot(
    d,
    aes(
      x = pos / 1e6,
      y = PC1,
      group = ind,
      color = POP
    )
  ) +
    geom_line(alpha = 0.4) +
    facet_wrap(~POP, ncol = 1,
               labeller = as_labeller(c(east = "ENP", west = "WNP"))) +
    scale_color_manual(values = POP_colors) +
    labs(
      title = paste0(cc),
      x = "Position (Mb)",
      y = "Window PC1"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      strip.text = element_text(size = 15, face = "bold"),
      legend.position = "none"
    )

  print(p)

  ggsave(
    filename = file.path(outdir, paste0(cc, "_sliding_PCA.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}
