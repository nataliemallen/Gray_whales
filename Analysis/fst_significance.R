# fst significance: block jackknife se + window subsampling (Table S2)
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")
set.seed(42)

w <- read.table("windows_AB.txt", col.names = c("chr","start","A","B","n"))
w <- w[w$B > 0, ]

weighted_fst <- function(A, B) sum(A) / sum(B)
global <- weighted_fst(w$A, w$B)

# ---- block jackknife (leave-one-window-out) ----
m  <- nrow(w)
jk <- sapply(seq_len(m), function(i) weighted_fst(w$A[-i], w$B[-i]))
jk_mean <- mean(jk)
jk_se   <- sqrt((m - 1) / m * sum((jk - jk_mean)^2))
cat(sprintf("Global weighted FST = %.5f   jackknife SE = %.2e\n", global, jk_se))

# block-bootstrap empirical p (resample 500 kb windows)
B <- 10000
boot <- replicate(B, { idx <- sample(m, m, replace = TRUE); weighted_fst(w$A[idx], w$B[idx]) })
p_boot <- (sum(boot <= 0) + 1) / (B + 1)
cat(sprintf("Block-bootstrap: min FST = %.5f over %d reps, empirical p = %.2g\n",
            min(boot), B, p_boot))

# or test Ho fst = 0
z    <- global / jk_se
logp <- log(2) + pnorm(-abs(z), log.p = TRUE)   # natural-log two-sided p
p    <- exp(logp)
cat(sprintf("Jackknife Z = %.1f   p = %.3g   (log10 p ~ %.0f)\n",
            z, p, logp / log(10)))
if (p < 2.2e-16) cat("  -> report as p < 2.2e-16 (or p < 0.001)\n")

# ---- window subsampling ----
sizes <- c(100, 500, 1000, 5000, 10000, 25000, 50000)
sizes <- sizes[sizes <= m]
res <- lapply(sizes, function(s) {
  vals <- replicate(1000, {
    idx <- sample(m, s)
    weighted_fst(w$A[idx], w$B[idx])
  })
  data.frame(n_windows = s, mean = mean(vals), sd = sd(vals),
             lo = quantile(vals, 0.025), hi = quantile(vals, 0.975))
})
res <- do.call(rbind, res)
print(res)
write.csv(res, "fst_subsampling_TableS2.csv", row.names = FALSE)
