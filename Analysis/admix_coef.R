# admixture-coefficient analyses: heterozygosity vs eastern ancestry (Fig S10, S11)
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")
library(ggplot2)

NSITES <- 3296954   # number of sites used to estimate admixture (from whale.beagle.gz)

# heterozygosity (from heterozygosity.R) and metadata
het  <- read.csv("heterozygosity_2026.csv", colClasses = c(ID = "character"))  # ID, POP, heterozygosity, ...
# K=2 admixture proportions, in bams_list.txt order
q    <- read.table("admix_K2_best.qopt")
meta <- read.table("sample_pop.tsv", header = TRUE, sep = "\t", colClasses = c(sample = "character"))
meta <- meta[meta$sample != "038046_ER-17-0237", ]
stopifnot(nrow(q) == nrow(meta))

# identify which column is "eastern" ancestry (higher mean within east)
east_idx <- which(meta$pop == "east")
east_col <- which.max(c(mean(q[east_idx,1]), mean(q[east_idx,2])))
admix <- data.frame(ID = meta$sample, POP = meta$pop, east_anc = q[, east_col])

# merge with H
m <- merge(admix, het[, c("ID","heterozygosity")], by = "ID")
wnp <- subset(m, POP == "west")

# 1b) Hump test: is heterozygosity ELEVATED at INTERMEDIATE ancestry?
# q(1-q) peaks at ancestry 0.5 and is the interpopulation-heterozygosity
# expectation under admixture; a positive, significant slope means heterozygosity
# is elevated in intermediate-ancestry individuals (i.e. recent admixture).
wnp$hybridity <- wnp$east_anc * (1 - wnp$east_anc)
hump <- lm(heterozygosity ~ hybridity, data = wnp)
cat("\nHump test  (heterozygosity ~ q(1-q)):\n"); print(summary(hump)$coefficients)

# equivalent quadratic form; a significant NEGATIVE quadratic term = a peak
quad <- lm(heterozygosity ~ east_anc + I(east_anc^2), data = wnp)
cat("\nQuadratic (heterozygosity ~ ancestry + ancestry^2):\n"); print(summary(quad)$coefficients)

# plot with the fitted hump (q(1-q)) curve
curve_df <- data.frame(east_anc = seq(min(wnp$east_anc), max(wnp$east_anc), length.out = 100))
curve_df$fit <- predict(hump, newdata = data.frame(hybridity = curve_df$east_anc * (1 - curve_df$east_anc)))
print(
  ggplot(wnp, aes(east_anc, heterozygosity)) +
    geom_line(data = curve_df, aes(east_anc, fit), color = "grey50") +
    geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey60") +
    geom_point(color = "#E7AF23", size = 2.5) +
    labs(x = "Eastern ancestry coefficient", y = "Heterozygosity",
         title = "WNP: heterozygosity vs admixture (grey = q(1-q) hump fit)") +
    theme_minimal()
)

# hump_pie.R -- WNP heterozygosity vs eastern ancestry, with each individual
# drawn as an admixture pie (eastern vs western K=2 proportions).
# Visual evidence that intermediate-ancestry (more evenly admixed) WNP
# individuals have higher heterozygosity, consistent with recent admixture.
# Needs: whale_k2_ordered.csv, heterozygosity_2026.csv
library(ggplot2)
library(scatterpie)          # install.packages("scatterpie")

# ---- eastern ancestry (K=2) ----
k2 <- read.csv("whale_k2_ordered.csv")
k2$id <- sub("\\.bam$", "", k2$SampleID)
east_is_c1 <- mean(k2$Cluster1[k2$POP == "east"]) > mean(k2$Cluster2[k2$POP == "east"])
k2$east_anc <- if (east_is_c1) k2$Cluster1 else k2$Cluster2

# ---- genome-wide heterozygosity, WNP only ----
het <- read.csv("heterozygosity_2026.csv", colClasses = c(ID = "character"))
w <- merge(k2[k2$POP == "west", c("id", "east_anc")],
           het[, c("ID", "heterozygosity")], by.x = "id", by.y = "ID")
w$eastern <- w$east_anc
w$western <- 1 - w$east_anc

# ---- fit the q(1-q) curve (peaks at intermediate ancestry) ----
fit <- lm(heterozygosity ~ I(east_anc * (1 - east_anc)), data = w)
cx  <- seq(0, 1, length.out = 200)
cyh <- predict(fit, data.frame(east_anc = cx))

# ---- map heterozygosity to a 0-1 axis so pies render as circles (coord_equal) ----
ylo <- min(w$heterozygosity); yhi <- max(w$heterozygosity)
pad <- 0.10 * (yhi - ylo); ylo <- ylo - pad; yhi <- yhi + pad
map_y  <- function(h) (h - ylo) / (yhi - ylo)
w$y    <- map_y(w$heterozygosity)
curve_df <- data.frame(x = cx, y = map_y(cyh))
ybreaks  <- pretty(c(min(w$heterozygosity), max(w$heterozygosity)), 4)

w$radius <- 0.018                 # pie size (in x units); tweak to taste
w$grp    <- seq_len(nrow(w))

p <- ggplot() +
  geom_line(data = curve_df, aes(x, y), color = "grey50") +
  geom_vline(xintercept = 0.5, linetype = "dotted", color = "grey60") +
  geom_scatterpie(aes(x = east_anc, y = y, group = grp, r = radius),
                  data = w, cols = c("eastern", "western"), color = NA) +
  scale_fill_manual(values = c(eastern = "#3FA5D5", western = "#E7AF23"),
                    name = "Ancestry") +
  scale_y_continuous(name = "Heterozygosity",
                     breaks = map_y(ybreaks), labels = ybreaks) +
  scale_x_continuous(name = "Eastern ancestry (NGSadmix K = 2)",
                     limits = c(-0.05, 1.05)) +
  coord_equal() +
  labs(title = "WNP heterozygosity is highest at intermediate ancestry") +
  theme_minimal()
print(p)
ggsave("hump_pie.png", p, width = 8, height = 6, dpi = 300)