# global dxy between ENP and WNP from per-population mafs (downsampled)
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")

e <- read.table(gzfile("east_maf.mafs.gz"), header = TRUE)
w <- read.table(gzfile("west_maf.mafs.gz"), header = TRUE)

key_e <- paste(e$chromo, e$position, sep = "_")
key_w <- paste(w$chromo, w$position, sep = "_")
shared <- intersect(key_e, key_w)

p1 <- e$knownEM[match(shared, key_e)]
p2 <- w$knownEM[match(shared, key_w)]

dxy_site <- p1 * (1 - p2) + p2 * (1 - p1)

# Denominator = TOTAL callable sites 
# number of sites analysed is the sum of the 2D SFS (joint SAF) from fst.sh.
nsites_total <- sum(scan("east_west.sfs", quiet = TRUE))

cat(sprintf("Polymorphic shared sites : %d\n", length(shared)))
cat(sprintf("Total callable sites     : %.0f\n", nsites_total))
cat(sprintf("Global Dxy = %.6f\n", sum(dxy_site, na.rm = TRUE) / nsites_total))
