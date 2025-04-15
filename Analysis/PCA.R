### R code for pca 
setwd("/Users/natal/R/whales")

install.packages("ggrepel")

library(ggrepel)
library(ggplot2)
library(dplyr)

# Load metadata dataframe
meta <- read.csv("whale_metadata_5-24.csv")  # Update with the actual path

# Load covariance matrix
turtles.cov <- as.matrix(read.table("whales.cov"))
turtles.cov <- as.matrix(read.table("LDpruned_whale_pca.cov"))

pca_result <- prcomp(turtles.cov, scale = TRUE)

# Combine PCA results with metadata
pca_df <- data.frame(pca_result$x)
pca_df <- cbind(meta, pca_df)  # Combining metadata and PCA results

# Extract % variation explained
var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

# Define the custom colors for east and west
POP_colors <- c("east" = "#1E88E5",  
                "west" = "#e51e8a")

# Create the PCA scatter plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) +
  geom_point() +
  labs(x = paste0("Principal Component 1 (", round(var_explained[1], 2), "%)"), 
       y = paste0("Principal Component 2 (", round(var_explained[2], 2), "%)"), 
       title = "PCA Scatter Plot") +
  
  # Apply the custom colors
  scale_color_manual(values = POP_colors) +
  
  # Customize theme for larger text, removing legend title, and light gray background
  theme_minimal(base_size = 16) +  # Set base size for axes and legend text
  theme(
    legend.title = element_blank(),  # Remove "POP" label from the legend
    legend.text = element_text(size = 14),  # Increase legend text size
    axis.title = element_text(size = 16),  # Increase axis title text size
    axis.text = element_text(size = 14),  # Increase axis label text size
    panel.background = element_rect(fill = "#f2f2f2", color = NA),  # Add light gray background
    #plot.background = element_rect(fill = "lightgray", color = NA)  # Match plot background with light gray
  )

# Print the plot
print(pca_plot)



#plot with file name labels instead
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = POP)) +
  geom_point() +
  geom_text_repel(aes(label = File), size = 3) +  # Add file name as a label
  labs(x = paste0("Principal Component 1 (", round(var_explained[1], 2), "%)"), 
       y = paste0("Principal Component 2 (", round(var_explained[2], 2), "%)"), 
       title = "PCA Scatter Plot")

print(pca_plot)

# Create scatter plot for PC2 vs PC3 with labels
pca_plot <- ggplot(pca_df, aes(x = PC2, y = PC3, color = POP)) +
  geom_point() +
  geom_text_repel(aes(label = File),  # Add labels from the "File" column
                  size = 3,           # Adjust label text size
                  box.padding = 0.35,  # Space around the text box
                  point.padding = 0.3) # Space between text and points

# Add axis labels with explained variance for PC2 and PC3
pca_plot <- pca_plot + labs(
  x = paste0("Principal Component 2 (", round(var_explained[2], 2), "%)"),
  y = paste0("Principal Component 3 (", round(var_explained[3], 2), "%)"),
  title = "PCA Scatter Plot (PC2 vs PC3)"
)

print(pca_plot)

# Create scatter plot for PC3 vs PC4 with labels
pca_plot <- ggplot(pca_df, aes(x = PC3, y = PC4, color = POP)) +
  geom_point() +
  geom_text_repel(aes(label = File),  # Add labels from the "File" column
                  size = 3,           # Adjust label text size
                  box.padding = 0.35,  # Space around the text box
                  point.padding = 0.3) # Space between text and points

# Add axis labels with explained variance for PC2 and PC3
pca_plot <- pca_plot + labs(
  x = paste0("Principal Component 3 (", round(var_explained[3], 2), "%)"),
  y = paste0("Principal Component 4 (", round(var_explained[4], 2), "%)"),
  title = "PCA Scatter Plot (PC2 vs PC3)"
)

print(pca_plot)

