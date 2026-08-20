# relatedness: r1/king degree classification and relatives-to-remove list (Fig 3, S1)
rm(list = ls())
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")

# NgsRelate v2 output (tab-separated, has a header with ida/idb/R1/KING)
kin <- read.table("whale_relate.tsv", header = TRUE, sep = "\t",
                  stringsAsFactors = FALSE, check.names = FALSE,
                  colClasses = c(ida = "character", idb = "character"))

# population assignment: columns sample_ID, pop  (east/west)
pop <- read.csv("whale_pop_2026.csv", header = TRUE, colClasses = c(sample_ID = "character"))   # cols: sample_ID, pop

# attach population labels for each member of the pair
kin$pop_a <- pop$pop[match(kin$ida, pop$sample_ID)]
kin$pop_b <- pop$pop[match(kin$idb, pop$sample_ID)]

# R1-KING category bounds (cols: R1_min, R1_max, KING_min, KING_max)
# 1st degree bound widened to include parent-offspring + full sibs (as in original)
bounds <- matrix(c(0.38,       Inf,       0.1767767, 0.3535534,   # 1st degree
                   0.20,       0.50,      0.08838835,0.1767767,    # 2nd degree
                   0.07142857, 0.4210526, 0.04419417,0.08838835,  # 3rd degree
                   -Inf,       Inf,      -Inf,       0.04419417),  # unrelated
                 ncol = 4, byrow = TRUE)
colors <- c("#E76BF3", "#00A5FF", "#00BC59", "#F8766D")
labs   <- c("1st degree", "2nd degree", "3rd degree", "Unrelated")

classify <- function(R1, KING) {
  for (i in 1:4) if (R1 >= bounds[i,1] & R1 <= bounds[i,2] &
                     KING >= bounds[i,3] & KING <= bounds[i,4]) return(labs[i])
  NA
}
kin$degree <- mapply(classify, kin$R1, kin$KING)
write.csv(kin, "kin_classified_2026.csv", row.names = FALSE)

plot_kin <- function(df, title) {
  
  pad.x <- diff(range(df$R1, na.rm = TRUE)) * 0.03
  pad.y <- diff(range(df$KING, na.rm = TRUE)) * 0.03
  
  plot(df$R1, df$KING,
       pch = 20, cex = 1, col = "black",
       xlim = range(df$R1, na.rm = TRUE) + c(-pad.x, pad.x),
       ylim = range(df$KING, na.rm = TRUE) + c(-pad.y, pad.y),
       xlab = "R1", ylab = "KING", main = title)
  
  for (i in 1:4) {
    sel <- df$R1 >= bounds[i,1] & df$R1 <= bounds[i,2] &
      df$KING >= bounds[i,3] & df$KING <= bounds[i,4]
    points(df$R1[sel], df$KING[sel],
           pch = 20, cex = 1, col = colors[i])
  }
  
  legend("bottomright", inset = c(0.02, 0.02),
         legend = labs,
         fill = colors,
         bty = "n",
         border = NA,
         pt.cex = 2,
         x.intersp = 2,
         y.intersp = 1.5)
}

plot_kin(kin, "All individuals")
plot_kin(kin[kin$pop_a == "east" & kin$pop_b == "east", ], "ENP")
plot_kin(kin[kin$pop_a == "west" & kin$pop_b == "west", ], "WNP")

kin_cross <- kin[
  (kin$pop_a == "east" & kin$pop_b == "west") |
    (kin$pop_a == "west" & kin$pop_b == "east"),
]

plot_kin(kin_cross, "East–West pairs")

tiff(
  filename = "kin_cross_EastWest.tif",
  width = 6,      # inches
  height = 5,     # inches
  units = "in",
  res = 300,      # 600 dpi (publication quality)
  compression = "lzw"
)

plot_kin(kin_cross, "East–West pairs")

dev.off()

remove <- c("038013_ER-16-0055",
            "038016_ER-16-0058")

# ---- relatives to remove for the relatives-excluded PCA ----
# Greedily remove one individual per related (<=3rd degree) pair to break all edges.
rel <- kin[kin$degree %in% c("1st degree","2nd degree","3rd degree"), c("ida","idb")]
remove <- character(0)
if (nrow(rel) > 0) {
  edges <- rel
  repeat {
    edges <- edges[!(edges$ida %in% remove | edges$idb %in% remove), ]
    if (nrow(edges) == 0) break
    tab <- sort(table(c(edges$ida, edges$idb)), decreasing = TRUE)
    remove <- c(remove, names(tab)[1])    # drop the most-connected individual
  }
}
writeLines(unique(remove), "relatives_remove.txt")
cat("Relatives removed for no-rel PCA:", length(unique(remove)), "\n")

# ---- verify field-identified cow-calf pairs ----
field_pairs <- data.frame(
  pair = 1:11,
  calf = c("ER_16_0054","ER_16_0061","ER_14_0160","ER_14_0173","ER_16_0069",
           "ER_16_0057","ER_16_0060","ER_16_0047","ER_16_0059","ER_16_0058","ER_16_0046"),
  mum  = c("ER_14_0147","ER_14_0151","ER_14_0159","ER_14_0159","ER_16_0045",
           "ER_16_0047","ER_16_0062","ER_16_0063","ER_16_0066","ER_16_0073","ER_14_0151"),
  stringsAsFactors = FALSE)

samples <- unique(c(kin$ida, kin$idb))
map_id  <- function(tok) {
  hit <- samples[grepl(gsub("_", "-", tok), samples, fixed = TRUE)]
  if (length(hit) == 1) hit else NA_character_       # NA if absent (e.g. removed) or ambiguous
}
field_pairs$calf_id <- vapply(field_pairs$calf, map_id, character(1))
field_pairs$mum_id  <- vapply(field_pairs$mum,  map_id, character(1))

# pull R1 / KING / classified degree for each pair (order-independent)
get_pair <- function(a, b) {
  if (is.na(a) || is.na(b)) return(data.frame(R1 = NA, KING = NA, degree = NA_character_))
  r <- kin[(kin$ida == a & kin$idb == b) | (kin$ida == b & kin$idb == a), ]
  if (nrow(r) == 0) return(data.frame(R1 = NA, KING = NA, degree = NA_character_))
  data.frame(R1 = r$R1[1], KING = r$KING[1], degree = r$degree[1])
}
pv <- do.call(rbind, Map(get_pair, field_pairs$calf_id, field_pairs$mum_id))
field_pairs <- cbind(field_pairs, pv)
field_pairs$verified_1st <- field_pairs$degree == "1st degree" & !is.na(field_pairs$degree)

write.csv(field_pairs, "cowcalf_verification.csv", row.names = FALSE)
cat(sprintf("Cow-calf pairs verified as first-degree: %d of %d\n",
            sum(field_pairs$verified_1st), nrow(field_pairs)))
not_found <- field_pairs[is.na(field_pairs$calf_id) | is.na(field_pairs$mum_id), "pair"]
if (length(not_found)) cat("Pairs with a missing member (not evaluable):", paste(not_found, collapse=", "), "\n")
print(field_pairs[, c("pair","calf","mum","R1","KING","degree","verified_1st")])

# ---- number of first-degree relative pairs by population ----
first  <- kin[kin$degree == "1st degree" & !is.na(kin$degree), ]
n_all  <- nrow(first)
n_east <- sum(first$pop_a == "east" & first$pop_b == "east")   # both ENP
n_west <- sum(first$pop_a == "west" & first$pop_b == "west")   # both WNP
n_btw  <- sum(first$pop_a != first$pop_b)                      # ENP x WNP
n_ind  <- function(df) length(unique(c(df$ida, df$idb)))       # unique individuals involved
cat(sprintf("First-degree pairs -- overall: %d | ENP: %d | WNP: %d | between: %d\n",
            n_all, n_east, n_west, n_btw))
cat(sprintf("Unique individuals in first-degree pairs -- overall: %d | ENP: %d | WNP: %d\n",
            n_ind(first),
            n_ind(first[first$pop_a=="east" & first$pop_b=="east", ]),
            n_ind(first[first$pop_a=="west" & first$pop_b=="west", ])))

