# vertical admixture bar plots for the tree figure (k=2-5)
setwd("/Users/natal/Documents/Purdue/Whale Project/2026_analysis")

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# vertical plots for tree 
# File names and corresponding K values
files <- c("whale_k2_ordered_tree.csv", "whale_k3_ordered_tree.csv", "whale_k4_ordered_tree.csv", "whale_k5_ordered_tree.csv")
k_values <- c(2, 3, 4, 5)

# Function to load and process data
load_and_process_data <- function(file, k) {
  # Load .csv dataframe
  qopt_data <- read.csv(file, header = TRUE)
  
  # Sort the dataframe based on the "Plot_Order" column
  qopt_data <- qopt_data[order(qopt_data$Plot_Order), ]
  
  # Convert cluster columns to numeric
  for (i in 1:k) {
    qopt_data[[paste0("Cluster", i)]] <- as.numeric(qopt_data[[paste0("Cluster", i)]])
  }
  
  # Reshape the data to long format
  plot_data <- qopt_data %>%
    pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster", values_to = "Proportion") %>%
    arrange(Plot_Order)  # Arrange the data by the "Plot_Order" column
  
  # Add the K value to the dataframe
  plot_data$k <- k
  
  return(plot_data)
}

# Load and process data for all K values
plot_data_list <- lapply(seq_along(files), function(i) load_and_process_data(files[i], k_values[i]))
plot_data <- bind_rows(plot_data_list)

# Define custom colors for the clusters
cluster_colors <- c("#3FA5D5","#E7AF23","#009e73","#d55e00","#713E5A")

# Create the vertical bar plot with facets arranged horizontally
admixture_plot <- ggplot(plot_data, aes(x = fct_rev(factor(Plot_Order)), y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_grid(cols = vars(k), scales = "free_x") + # Arrange K values horizontally
  coord_flip() +  # Rotate the bars
  theme(axis.text.y = element_blank(),  # Remove y-axis text (previous x-axis)
        axis.ticks.y = element_blank(), # Remove y-axis ticks (previous x-axis)
        strip.background = element_blank(),  # Clean up facet labels
        strip.text = element_text(size = 12)) + # Adjust facet label text size
  scale_fill_manual(values = cluster_colors) +
  labs(x = "Samples", y = "Proportion")  # Adjust axis labels

# Display the plot
print(admixture_plot)


