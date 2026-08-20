# genome-wide heterozygosity: ENP vs WNP and cross-species (Fig 6)
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")
library(ggplot2); library(dplyr)
POP_colors <- c("east" = "#3FA5D5", "west" = "#E7AF23")

# metadata: sample, pop, breadth, depth  (matches lists/sample_pop.tsv, 038046 removed)
meta <- read.table("sample_pop.tsv", header = TRUE, sep = "\t", colClasses = c(sample = "character"))
meta <- meta[meta$sample != "038046_ER-17-0237", ]

# read each *_est.ml; folded single-sample SFS -> H = het / total
het <- sapply(meta$sample, function(s) {
  f <- paste0("het/", s, "_est.ml")
  v <- scan(f, quiet = TRUE)
  v[2] / sum(v)
})
data <- data.frame(ID = meta$sample, POP = meta$pop,
                   heterozygosity = as.numeric(het),
                   breadth = meta$breadth, depth = meta$depth)
write.csv(data, "heterozygosity_2026.csv", row.names = FALSE)
cat(sprintf("Mean genome-wide H = %.5f\n", mean(data$heterozygosity)))

cat(sprintf("Mean genome-wide H = %.5f\n", mean(data$heterozygosity)))

data %>%
  group_by(POP) %>%
  summarise(
    n = n(),
    mean_H = mean(heterozygosity),
    sd_H = sd(heterozygosity)
  ) %>%
  print()

# boxplot by site (original style)
data$POP <- factor(data$POP, levels = c("east", "west"))

p <- ggplot(data, aes(POP, heterozygosity, color = POP)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 3) +
  scale_color_manual(
    values = POP_colors,
    labels = c("ENP", "WNP")
  ) +
  scale_x_discrete(labels = c("east" = "ENP", "west" = "WNP")) +
  labs(
    x = "Sample Site",
    y = "Heterozygosity",
    color = "Sample Site"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

print(p)

# Load the necessary libraries
library(ggplot2)
library(dplyr)

# Create a data frame with the heterozygosity values and species labels
data <- data.frame(
  Species = factor(c("Blue whale", "Fin whale", "Humpback whale", 
                     "Minke whale", "N.At. right whale", "Sei whale A", "Sei whale B", 
                     "Árnason Gray whale", "ENP Gray whale", "WNP Gray whale"),
                   levels = c("Blue whale", "Fin whale", "Humpback whale", 
                              "Minke whale", "N.At. right whale", "Sei whale A", "Sei whale B", 
                              "Árnason Gray whale", "ENP Gray whale", "WNP Gray whale")),
  heterozygosity = c(0.00295, 0.00185, 0.0012, 
                     0.0006, 0.00175, 0.00075, 0.00077, 
                     0.0007, 0.000460, 0.000483)
)

data <- data %>%
  mutate(Species = reorder(Species, -heterozygosity))

# Create the bar plot with different colors for each species
ggplot(data, aes(x = Species, y = heterozygosity, fill = Species)) +
  geom_bar(stat = "identity") +
  theme_minimal(base_size = 16) +  # increases all text globally
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title  = element_text(size = 16),
    plot.title  = element_text(size = 18, face = "bold"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Species", y = "Heterozygosity") +
  scale_fill_manual(values = c("Blue whale" = "skyblue", 
                               "Fin whale" = "#316d16", 
                               "Árnason Gray whale" = "#cc79a7", 
                               "Humpback whale" = "#d55e00", 
                               "Minke whale" = "#b87333", 
                               "N.At. right whale" = "#9882cf", 
                               "Sei whale A" = "turquoise", 
                               "Sei whale B" = "#f0e442", 
                               "ENP Gray whale" = "#3FA5D5", 
                               "WNP Gray whale" = "#E7AF23"))



# Wilcoxon ENP vs WNP
print(wilcox.test(heterozygosity ~ POP, data = data))

# regression on breadth and depth
cat("H ~ breadth:\n"); print(summary(lm(heterozygosity ~ breadth, data = data))$coefficients)
cat("H ~ depth:\n");   print(summary(lm(heterozygosity ~ depth,   data = data))$coefficients)

# ---- cross-species comparison (Fig 6B) ----
# species_het.csv: species, heterozygosity  (gray whale + other rorquals/marine mammals)
if (file.exists("species_het.csv")) {
  sp <- read.csv("species_het.csv")
  sp <- sp[order(sp$heterozygosity), ]
  sp$species <- factor(sp$species, levels = sp$species)
  print(
    ggplot(sp, aes(species, heterozygosity)) +
      geom_col(fill = "grey40") +
      coord_flip() +
      labs(x = NULL, y = "Heterozygosity") + theme_minimal()
  )
}
