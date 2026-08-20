# hybrid-index vs interclass-heterozygosity triangle for f1 detection
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")
library(ggplot2)

# ---- inputs ----
aims <- read.table("aims.tsv", header = TRUE, sep = "\t", colClasses = c(chr = "character"))
inds <- readLines("tri_geno.inds")
pop  <- read.csv("whale_pop_2026.csv", colClasses = c(sample_ID = "character"))  # sample_ID, pop
ref6 <- sub("\\.bam$", "", basename(readLines("wnp6_extreme_bams.txt")))         # western reference IDs

# genotype posteriors: cols = chr, pos, then (P0,P1,P2) per individual
g <- read.table(gzfile("tri_geno.geno.gz"), header = FALSE, stringsAsFactors = FALSE)
g[[1]] <- as.character(g[[1]])

# align posterior rows to the AIM table (by chr:pos) to get polarization
idx  <- match(paste(g[[1]], g[[2]]), paste(aims$chr, aims$pos))
keep <- !is.na(idx)
g <- g[keep, ]; aim <- aims[idx[keep], ]
east_is_minor <- aim$enp_freq > aim$wnp_freq     # is the eastern allele the minor allele?

# ---- per-individual hybrid index (HI) and interclass heterozygosity (IH) ----
nInd <- length(inds)
HI <- IH <- nSNP <- numeric(nInd)
for (i in seq_len(nInd)) {
  b  <- 2 + (i - 1) * 3
  P0 <- g[[b + 1]]; P1 <- g[[b + 2]]; P2 <- g[[b + 3]]
  # drop no-data sites: uniform posterior (~1/3 each) under the uniform prior
  informative <- !(abs(P0 - 1/3) < 1e-3 & abs(P1 - 1/3) < 1e-3 & abs(P2 - 1/3) < 1e-3)
  east_frac <- ifelse(east_is_minor, (P1 + 2 * P2) / 2, (P1 + 2 * P0) / 2)
  HI[i]   <- mean(east_frac[informative], na.rm = TRUE)
  IH[i]   <- mean(P1[informative],        na.rm = TRUE)
  nSNP[i] <- sum(informative, na.rm = TRUE)
}
df <- data.frame(ind = inds, HI = HI, IH = IH, nSNP = nSNP,
                 pop = pop$pop[match(inds, pop$sample_ID)])
df$group <- ifelse(df$pop == "east", "ENP",
             ifelse(df$ind %in% ref6, "WNP (reference)", "WNP (test)"))
write.csv(df, "triangle_HI_IH.csv", row.names = FALSE)
cat(sprintf("AIMs used: %d   individuals: %d\n", nrow(aim), nInd))

# ---- plot ----
tri <- data.frame(x = c(0, 1, 0.5, 0), y = c(0, 0, 1, 0))          # triangle outline
cls <- data.frame(class = c("ENP", "WNP", "F1", "F2", "BC-ENP", "BC-WNP"),
                  HI = c(1, 0, 0.5, 0.5, 0.75, 0.25),
                  IH = c(0, 0, 1.0, 0.5, 0.50, 0.50))
cols <- c("ENP" = "#3FA5D5", "WNP (reference)" = "#8a6200", "WNP (test)" = "#E7AF23")

p <- ggplot() +
  geom_path(data = tri, aes(x, y), linetype = "dashed", color = "grey60") +
  geom_point(data = cls, aes(HI, IH), shape = 4, size = 3, stroke = 1) +
  geom_text(data = cls, aes(HI, IH, label = class), vjust = -0.9, size = 3) +
  geom_point(data = df, aes(HI, IH, color = group), size = 2.6, alpha = 0.85) +
  scale_color_manual(values = cols, name = NULL) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.02)) +
  labs(x = "Hybrid index (eastern ancestry)",
       y = "Interclass heterozygosity",
       title = "Hybrid index vs interclass heterozygosity") +
  theme_minimal()
print(p)
ggsave("triangle_plot.png", p, width = 7, height = 6, dpi = 300)

# ---- NewHybrids input (hard calls, posterior >= 0.8; 0 = missing) ----
# Two-allele coding: 11 = hom major, 12 = het, 22 = hom minor.
call <- matrix(0L, nrow = nInd, ncol = nrow(aim))
for (i in seq_len(nInd)) {
  b  <- 2 + (i - 1) * 3
  P  <- cbind(g[[b + 1]], g[[b + 2]], g[[b + 3]])
  mx <- max.col(P, ties.method = "first"); conf <- P[cbind(seq_len(nrow(P)), mx)]
  code <- c(11L, 12L, 22L)[mx]; code[conf < 0.8] <- 0L
  call[i, ] <- code
}
con <- file("newhybrids_input.txt", "w")
writeLines(c(sprintf("NumIndivs %d", nInd),
             sprintf("NumLoci %d",  nrow(aim)),
             "Digits 1", "Format Lumped",
             paste("LocusNames", paste0("L", seq_len(nrow(aim)), collapse = " "))), con)
for (i in seq_len(nInd))
  writeLines(paste(inds[i], paste(call[i, ], collapse = " ")), con)
close(con)
cat("Wrote triangle_plot.png, triangle_HI_IH.csv, newhybrids_input.txt\n")

# replot for an interpretable figure

d    <- read.csv("triangle_HI_IH.csv")
meta <- read.table("sample_pop.tsv", header = TRUE, sep = "\t",
                   colClasses = c(sample = "character"))
d$depth <- meta$depth[match(d$ind, meta$sample)]

# ---- rescale ancestry axis to empirical anchors (0 = western ref, 1 = ENP) ----
east_anchor <- mean(d$HI[d$group == "ENP"])
west_anchor <- mean(d$HI[d$group == "WNP (reference)"])
rescale <- function(hi) (hi - west_anchor) / (east_anchor - west_anchor)
d$HIr <- rescale(d$HI)

# ---- depth confound: regress IH on depth, keep residual (depth-corrected IH) ----
fit <- lm(IH ~ depth, data = d)
d$IH_resid <- resid(fit)
cat(sprintf("IH ~ depth: slope p = %.2g, R2 = %.2f  (Spearman shown separately)\n",
            summary(fit)$coefficients[2, 4], summary(fit)$r.squared))

# ---- class expectations calibrated to the ACTUAL AIM allele frequencies ----
a   <- read.table("aims.tsv", header = TRUE, sep = "\t", colClasses = c(chr = "character"))
eim <- a$enp_freq > a$wnp_freq
pe  <- ifelse(eim, a$enp_freq, 1 - a$enp_freq)     # eastern-allele freq in ENP
pw  <- ifelse(eim, a$wnp_freq, 1 - a$wnp_freq)     # eastern-allele freq in WNP
q   <- (pe + pw) / 2
cls <- data.frame(
  class = c("ENP","WNP","F1","F2","BC->ENP","BC->WNP"),
  HI = c(mean(pe), mean(pw), mean(q), mean(q), mean((q+pe)/2), mean((q+pw)/2)),
  IH = c(mean(2*pe*(1-pe)), mean(2*pw*(1-pw)),
         mean(pe*(1-pw)+(1-pe)*pw), mean(2*q*(1-q)),
         mean(q*(1-pe)+(1-q)*pe), mean(q*(1-pw)+(1-q)*pw)))
# calibrate theoretical IH to the observed het scale using the pure classes
k <- mean(d$IH[d$group %in% c("ENP","WNP (reference)")]) /
  mean(cls$IH[cls$class %in% c("ENP","WNP")])
cls$IHc <- cls$IH * k
cls$HIr <- rescale(cls$HI)
cat(sprintf("Het calibration factor k = %.2f\n", k)); print(cls)

cols <- c("ENP" = "#3FA5D5", "WNP (reference)" = "#8a6200", "WNP (test)" = "#E7AF23")

# ---- Plot 1: calibrated triangle, point size = depth (makes the confound visible) ----
p1 <- ggplot() +
  geom_point(data = d, aes(HIr, IH, color = group, size = depth), alpha = 0.8) +
  geom_point(data = cls, aes(HIr, IHc), shape = 4, size = 3, stroke = 1.1) +
  geom_text(data = cls, aes(HIr, IHc, label = class), vjust = -0.9, size = 3) +
  scale_color_manual(values = cols, name = NULL) +
  scale_size_continuous(range = c(1, 5), name = "depth (X)") +
  labs(x = "Relative eastern ancestry (0 = western ref, 1 = ENP)",
       y = "Interclass heterozygosity",
       title = "Calibrated triangle (marker size = sequencing depth)") +
  theme_minimal()
print(p1); ggsave("triangle_calibrated.png", p1, width = 8, height = 6, dpi = 300)

# ---- Plot 2: depth-corrected hybrid test (the honest one) ----
# A recent F1 should sit near x = 0.5 AND clearly above 0 on the y-axis
# (excess heterozygosity beyond what depth explains).
p2 <- ggplot(d, aes(HIr, IH_resid, color = group)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey60") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey60") +
  geom_point(size = 2.6, alpha = 0.85) +
  scale_color_manual(values = cols, name = NULL) +
  labs(x = "Relative eastern ancestry (0 = western ref, 1 = ENP)",
       y = "Interclass heterozygosity (depth-corrected residual)",
       title = "Depth-corrected: excess heterozygosity vs ancestry") +
  theme_minimal()
print(p2); ggsave("triangle_depthcorrected.png", p2, width = 7.5, height = 6, dpi = 300)

# ---- rank candidates: elevated depth-corrected het, near mid-ancestry ----
cat("\nTop individuals by depth-corrected interclass heterozygosity:\n")
print(head(d[order(-d$IH_resid), c("ind","group","depth","HIr","IH","IH_resid")], 8))

# hybrid_detection.R -- FINAL hybrid-detection figure (replaces the binomial-band method)
# hybrid_detection.R -- interclass heterozygosity vs eastern ancestry
# x = eastern ancestry (NGSadmix K = 2): a model-based ancestry estimate that,
#     unlike a dosage-based hybrid index, is not inflated by heterozygosity.
# y = interclass heterozygosity at ancestry-informative markers, depth-adjusted.
# Grey line +/- 2 SD band = expected heterozygosity as a function of ancestry,
# fit on the (pure) reference individuals. A recent F1 would sit near ancestry
# 0.5 with heterozygosity well ABOVE this expectation.

library(ggplot2)

# eastern ancestry from admixture
k2 <- read.csv("whale_k2_ordered.csv")
k2$id <- sub("\\.bam$", "", k2$SampleID)
east_is_c1 <- mean(k2$Cluster1[k2$POP == "east"]) > mean(k2$Cluster2[k2$POP == "east"])
k2$east_anc <- if (east_is_c1) k2$Cluster1 else k2$Cluster2

# interclass heterozygosity + depth adjustment
tri  <- read.csv("triangle_HI_IH.csv")
meta <- read.table("sample_pop.tsv", header = TRUE, sep = "\t", colClasses = c(sample = "character"))
tri$depth    <- meta$depth[match(tri$ind, meta$sample)]
tri$IH_adj   <- resid(lm(IH ~ depth, data = tri)) + mean(tri$IH)     # remove depth confound
tri$east_anc <- k2$east_anc[match(tri$ind, k2$id)]
tri$group <- ifelse(tri$group == "ENP", "ENP",
                    ifelse(tri$group == "WNP (reference)", "WNP (low-eastern ref.)", "WNP"))
cols <- c("ENP" = "#3FA5D5", "WNP" = "#E7AF23", "WNP (low-eastern ref.)" = "#8a6200")

# expected heterozygosity vs ancestry, fit on PURE (reference) individuals only
pure <- tri$group != "WNP"
fit  <- lm(IH_adj ~ east_anc, data = tri[pure, ])
sdp  <- sd(resid(fit))
band <- data.frame(east_anc = seq(0, 1, 0.01))
band$exp <- predict(fit, band)

# flag individuals >2 SD above the ancestry expectation
tri$excess_z <- (tri$IH_adj - predict(fit, tri)) / sdp

p <- ggplot() +
  geom_ribbon(data = band, aes(east_anc, ymin = exp - 2*sdp, ymax = exp + 2*sdp),
              fill = "grey85", alpha = 0.6) +
  geom_line(data = band, aes(east_anc, exp), color = "grey50") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey60") +
  geom_point(data = tri, aes(east_anc, IH_adj, color = group), size = 2.6, alpha = 0.85) +
  scale_color_manual(values = cols, name = NULL) +
  scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Eastern ancestry (NGSadmix K = 2)",
       y = "Interclass heterozygosity (depth-adjusted)",
       title = "Interclass heterozygosity vs ancestry",
       subtitle = "grey band = pure-individual expectation ± 2 SD; dotted line = F1 ancestry (0.5)") +
  theme_minimal()
print(p)
ggsave("hybrid_detection.png", p, width = 7.5, height = 6, dpi = 300)

# report individuals above the expectation (potential recent-admixture signal)
cat("Individuals >2 SD above the ancestry expectation:\n")
print(tri[order(-tri$excess_z), c("ind","group","east_anc","IH_adj","excess_z")][1:6, ])