### R code for pca 
setwd("/Users/natal/R/whales")

library(ggplot2)
library(dplyr)

# Load metadata dataframe
meta <- read.csv("whale_metadata_5-24.csv")  # Update with the actual path

# Load covariance matrix
turtles.cov <- as.matrix(read.table("whales.cov"))

pca_result <- prcomp(turtles.cov, scale = TRUE)

# Combine PCA results with metadata
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)  # Combining metadata and PCA results

# Extract % variation explained
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Create scatter plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) +
  geom_point() +
  labs(x = paste0("Principal Component 1 (", round(var_explained[1], 2), "%)"), 
       y = paste0("Principal Component 2 (", round(var_explained[2], 2), "%)"), 
       title = "PCA Scatter Plot")

print(pca_plot)
