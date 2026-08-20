# folded sfs and tajima's d stats: ENP vs WNP (Fig 5C, S10)
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")
set.seed(42)
library(ggplot2); library(dplyr); library(tidyr)
POP_colors <- c("east" = "#3FA5D5", "west" = "#E7AF23")

# ---------- folded SFS (Fig 5C) ----------
# realSFS writes a full-length (2N+1) vector; for a FOLDED spectrum only the
# minor-allele classes 0..N carry counts and the upper half are zeros. Keep just
# the informative classes 1..N (drop the monomorphic class 0 and the folded zeros).
read_sfs <- function(f) {
  s <- scan(f, quiet = TRUE)
  N <- (length(s) - 1) / 2      # number of diploid individuals -> N folded classes
  s[2:(N + 1)]                  # minor-allele-count classes 1..N (raw counts)
}
e_sfs <- read_sfs("east.sfs"); w_sfs <- read_sfs("west.sfs")
sfs_df <- data.frame(freq = seq_along(e_sfs),
                     east = e_sfs, west = w_sfs) %>%
  pivot_longer(c(east, west), names_to = "POP", values_to = "nsites")
print(
  ggplot(sfs_df, aes(freq, nsites, fill = POP)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = POP_colors) +
    labs(x = "Minor allele count", y = "Number of sites",
         title = "Folded SFS") + theme_minimal() +
    theme(
      axis.title = element_text(size = 16)
    )
)

# KS test: normalise to proportions internally, then expand into pseudo-observations
expand_sfs <- function(x, n = 1e5) { p <- x / sum(x); rep(seq_along(p), round(p * n)) }
ks <- ks.test(expand_sfs(e_sfs), expand_sfs(w_sfs))
cat(sprintf("KS test: D = %.4f, p = %.3g\n", ks$statistic, ks$p.value))

# ---------- Tajima's D (Fig S10) ----------
# thetaStat .pestPG columns include Chr, WinCenter, ..., Tajima, ..., nSites
read_pestpg <- function(f) {
  d <- read.table(f, header = TRUE, comment.char = "")
  names(d) <- gsub("^X\\.+", "", names(d))
  d
}
e_taj <- read_pestpg("east_50kb.pestPG")$Tajima
w_taj <- read_pestpg("west_50kb.pestPG")$Tajima
e_taj <- e_taj[is.finite(e_taj)]; w_taj <- w_taj[is.finite(w_taj)]

wt <- wilcox.test(e_taj, w_taj)
cat(sprintf("Tajima's D  median ENP = %.3f  WNP = %.3f  (diff = %.3f)\n",
            median(e_taj), median(w_taj), median(e_taj) - median(w_taj)))
cat(sprintf("Wilcoxon rank-sum p = %.3g\n", wt$p.value))

boot_ci <- function(x, R = 10000) {
  m <- replicate(R, mean(sample(x, replace = TRUE)))
  quantile(m, c(0.025, 0.975))
}
cat("Bootstrap 95% CI mean Tajima's D  ENP:", round(boot_ci(e_taj), 4), "\n")
cat("Bootstrap 95% CI mean Tajima's D  WNP:", round(boot_ci(w_taj), 4), "\n")

taj_df <- rbind(data.frame(Tajima = e_taj, POP = "east"),
                data.frame(Tajima = w_taj, POP = "west"))
# print(
#   ggplot(taj_df, aes(POP, Tajima, color = POP)) +
#     geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.2) +
#     scale_color_manual(values = POP_colors) +
#     labs(x = "Sample Site", y = "Tajima's D (50 kb windows)") + theme_minimal()
# )

taj_df <- rbind(
  data.frame(Tajima = e_taj, POP = "ENP"),
  data.frame(Tajima = w_taj, POP = "WNP")
)

taj_df$POP <- factor(taj_df$POP, levels = c("ENP", "WNP"))

print(
  ggplot(taj_df, aes(x = Tajima, fill = POP, color = POP)) +
    geom_density(alpha = 0.4, linewidth = 1) +
    scale_fill_manual(
      values = c(
        "ENP" = "#3FA5D5",
        "WNP" = "#E7AF23"
      )
    ) +
    scale_color_manual(
      values = c(
        "ENP" = "#3FA5D5",
        "WNP" = "#E7AF23"
      )
    ) +
    labs(
      x = "Tajima's D (50 kb windows)",
      y = "Density",
      fill = "Sample Site",
      color = "Sample Site"
    ) +
    scale_x_continuous(limits = c(-3, 3)) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )
)

# updated
# SFS + Tajima's D statistics (Figs 5C, S10)
# - KS test on folded SFS distributions
# - Wilcoxon rank-sum on windowed Tajima's D + median difference (effect size)
# - bootstrap 95% CIs (10,000 reps) for mean Tajima's D per population
set.seed(42)
library(ggplot2); library(dplyr); library(tidyr)
POP_colors <- c("east" = "#3FA5D5", "west" = "#E7AF23")

# ---------- folded SFS (Fig 5C) ----------
# realSFS writes a full-length (2N+1) vector; for a FOLDED spectrum only the
# minor-allele classes 0..N carry counts and the upper half are zeros. Keep just
# the informative classes 1..N (drop the monomorphic class 0 and the folded zeros).
read_sfs <- function(f) {
  s <- scan(f, quiet = TRUE)
  N <- (length(s) - 1) / 2      # number of diploid individuals -> N folded classes
  s[2:(N + 1)]                  # minor-allele-count classes 1..N (raw counts)
}
e_sfs <- read_sfs("east.sfs"); w_sfs <- read_sfs("west.sfs")
sfs_df <- data.frame(freq = seq_along(e_sfs),
                     east = e_sfs, west = w_sfs) %>%
  pivot_longer(c(east, west), names_to = "POP", values_to = "nsites")
print(
  ggplot(sfs_df, aes(freq, nsites, fill = POP)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = POP_colors) +
    labs(x = "Minor allele frequency class", y = "Number of sites",
         title = "Folded SFS") + theme_minimal()
)

# SFS divergence, reported as an EFFECT SIZE (Kolmogorov-Smirnov D statistic =
# the maximum difference between the two cumulative spectra). We do NOT attach a
# p-value: the earlier approach fabricated pseudo-observations, which makes the KS
# p-value scale with an arbitrary n; and with genome-wide site counts any nonzero
# difference is "significant", so the magnitude of D (0 = identical) is the
# informative quantity.
e_prop <- e_sfs / sum(e_sfs)
w_prop <- w_sfs / sum(w_sfs)
ks_D   <- max(abs(cumsum(e_prop) - cumsum(w_prop)))
cat(sprintf("SFS divergence: KS D = %.4f  (0 = identical distributions)\n", ks_D))

# Optional legitimate test on the REAL per-class site counts (G-test of
# independence). Reported for completeness only -- with millions of sites it is
# essentially always significant, so interpret D (above), not this p.
gtest_tab <- rbind(east = e_sfs, west = w_sfs)
sup <- suppressWarnings(chisq.test(gtest_tab))
cat(sprintf("  (reference chi-square on class counts: X2 = %.1f, p = %.3g; large-n, see D)\n",
            sup$statistic, sup$p.value))

# ---------- Tajima's D (Fig S10) ----------
# thetaStat .pestPG columns include Chr, WinCenter, ..., Tajima, ..., nSites
read_pestpg <- function(f) {
  d <- read.table(f, header = TRUE, comment.char = "")
  names(d) <- gsub("^X\\.+", "", names(d))
  d
}
e_taj <- read_pestpg("east_50kb.pestPG")$Tajima
w_taj <- read_pestpg("west_50kb.pestPG")$Tajima
e_taj <- e_taj[is.finite(e_taj)]; w_taj <- w_taj[is.finite(w_taj)]

wt <- wilcox.test(e_taj, w_taj)
cat(sprintf("Tajima's D  median ENP = %.3f  WNP = %.3f  (diff = %.3f)\n",
            median(e_taj), median(w_taj), median(e_taj) - median(w_taj)))
cat(sprintf("Wilcoxon rank-sum p = %.3g\n", wt$p.value))

boot_ci <- function(x, R = 10000) {
  m <- replicate(R, mean(sample(x, replace = TRUE)))
  quantile(m, c(0.025, 0.975))
}
cat("Bootstrap 95% CI mean Tajima's D  ENP:", round(boot_ci(e_taj), 4), "\n")
cat("Bootstrap 95% CI mean Tajima's D  WNP:", round(boot_ci(w_taj), 4), "\n")

taj_df <- rbind(data.frame(Tajima = e_taj, POP = "east"),
                data.frame(Tajima = w_taj, POP = "west"))
print(
  ggplot(taj_df, aes(POP, Tajima, color = POP)) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, alpha = 0.2) +
    scale_color_manual(values = POP_colors) +
    labs(x = "Sample Site", y = "Tajima's D (50 kb windows)") + theme_minimal()
)