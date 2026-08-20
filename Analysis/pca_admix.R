# pca with k=2 admixture pie charts

library(scatterpie)

setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# File names and corresponding K values
files <- c("whale_k2_ordered.csv", "whale_k3_ordered.csv", "whale_k4_ordered.csv", "whale_k5_ordered.csv")
k_values <- c(2, 3, 4, 5)

# Function to load and process data
load_and_process_data <- function(file, k) {
  # Load .csv dataframe
  qopt_data <- read.csv(file, header = TRUE)
  
  # Sort the dataframe based on the "Order" column
  qopt_data <- qopt_data[order(qopt_data$Order), ]
  
  # Convert cluster columns to numeric
  for (i in 1:k) {
    qopt_data[[paste0("Cluster", i)]] <- as.numeric(qopt_data[[paste0("Cluster", i)]])
  }
  
  # Reshape the data to long format
  plot_data <- qopt_data %>%
    pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster", values_to = "Proportion") %>%
    arrange(Order)  # Arrange the data by the "Order" column
  
  # Add the K value to the dataframe
  plot_data$k <- k
  
  return(plot_data)
}

# Load and process data for all K values
plot_data_list <- lapply(seq_along(files), function(i) load_and_process_data(files[i], k_values[i]))
plot_data <- bind_rows(plot_data_list)

# Filter for K=2
admix_k2 <- plot_data %>%
  filter(k == 2) %>%
  select(Individual = Order, Cluster, Proportion) %>%
  pivot_wider(names_from = Cluster, values_from = Proportion)

library(ggrepel)
library(ggplot2)
library(dplyr)

# Load metadata dataframe
meta <- read.csv("whale_metadata_2026.csv")  # Update with the actual path

# Load covariance matrix
whales.cov <- as.matrix(read.table("whale.cov"))
#whales.cov <- as.matrix(read.table("LDpruned_whale_pca.cov"))

pca_result <- prcomp(whales.cov, scale = TRUE)

# Compute % variance explained by each PC
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

# Combine PCA results with metadata
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)

# Merge PCA coordinates with K=2 admixture data
pca_admix <- pca_df %>%
  left_join(admix_k2, by = c("Order" = "Individual"))

pca_admix

# PCA with pie charts
pca_pie <- ggplot() +
  geom_scatterpie(data = pca_admix, 
                  aes(x = PC1, y = PC2, r = 0.3),  # r = radius of pie
                  cols = c("Cluster1", "Cluster2")) +
  coord_equal() +
  labs(x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(var_explained[2], 2), "%)"),
       title = "PCA with K=2 Admixture Proportions") +
  scale_fill_manual(values = c("Cluster1" = "#3FA5D5", "Cluster2" = "#E7AF23")) +
  theme_minimal(base_size = 12)

print(pca_pie)

#
library(scatterpie)

# Merge PCA coordinates with K=2 admixture data
pca_admix <- pca_df %>%
  left_join(admix_k2, by = c("Order" = "Individual"))

library(scatterpie)

# Merge PCA coordinates with K=2 admixture data
pca_admix <- pca_df %>%
  left_join(admix_k2, by = c("Order" = "Individual"))

# Define colors
cluster_colors <- c("Cluster1" = "#3FA5D5", "Cluster2" = "#E7AF23")
POP_colors <- c("east" = "#A06CD5", "west" = "#EDE5A6")  # bright, contrasting

# PCA with pies + population-colored borders
pca_pie <- ggplot() +
  geom_scatterpie(
    data = pca_admix,
    aes(x = PC1, y = PC2, r = 0.15, group = Order, colour = POP), # <- FIXED here
    cols = c("Cluster1", "Cluster2"),
    size = 1
  ) +
  coord_equal() +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)"),
    title = "PCA with K=2 Admixture Proportions"
  ) +
  scale_fill_manual(values = cluster_colors, name = "Admixture") +
  scale_color_manual(values = POP_colors, name = "Population") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "right")

print(pca_pie)

library(dplyr)

# Recode POP for nicer facet labels
pca_admix <- pca_admix %>%
  mutate(POP_label = recode(POP,
                            "east" = "ENP",
                            "west" = "WNP"))

pca_admix <- pca_admix %>%
  mutate(
    POP_label = factor(POP_label, levels = rev(c("ENP", "WNP")))
  )


# Then use POP_label in the facet
pca_pie <- ggplot() +
  geom_scatterpie(
    data = pca_admix,
    aes(x = PC1, y = PC2, r = 0.2, group = Order),
    cols = c("Cluster1", "Cluster2"),
    color = "black", size = 0.3
  ) +
  coord_equal() +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)"),
    title = "PCA with K=2 Admixture Proportions"
  ) +
  facet_wrap(~POP_label, ncol = 1, scales = "fixed") +   # use new labels
  scale_fill_manual(values = cluster_colors, name = "Admixture") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "right")

print(pca_pie)

ggsave(filename = "pca_admix.jpg", plot=pca_pie, dpi = 300)


# East-only
pca_pie_east <- ggplot() +
  geom_scatterpie(
    data = filter(pca_admix, POP == "east"),
    aes(x = PC1, y = PC2, r = 0.06, group = Order),
    cols = c("Cluster1", "Cluster2"),
    color = "black", size = 0.5
  ) +
  coord_equal() +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)"),
    title = "PCA with K=2 Admixture Proportions (East)"
  ) +
  scale_fill_manual(values = cluster_colors, name = "Admixture") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "right")

print(pca_pie_east)

# West-only
pca_pie_west <- ggplot() +
  geom_scatterpie(
    data = filter(pca_admix, POP == "west"),
    aes(x = PC1, y = PC2, r = 0.15, group = Order),
    cols = c("Cluster1", "Cluster2"),
    color = "black", size = 0.5
  ) +
  coord_equal() +
  labs(
    x = paste0("PC1 (", round(var_explained[1], 2), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 2), "%)"),
    title = "PCA with K=2 Admixture Proportions (West)"
  ) +
  scale_fill_manual(values = cluster_colors, name = "Admixture") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "right")

print(pca_pie_east)
print(pca_pie_west)

