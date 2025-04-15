### R code for pca 
setwd("/Users/natal/R/whales")

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Function to load and plot data for a specific K value
plot_admixture <- function(K) {
  # Load .csv dataframe
  file_name <- paste0("whale_k", K, "_combined.csv")
  qopt_data <- read.csv(file_name, header = TRUE)
  
  # Convert columns to numeric if needed
  for (i in 1:K) {
    qopt_data[[paste0("Cluster", i)]] <- as.numeric(qopt_data[[paste0("Cluster", i)]])
  }
  
  # Sort the dataframe based on the "Order" column
  qopt_data <- qopt_data[order(qopt_data$Order), ]
  
  # Create a sequence from 1 to n (number of unique individuals)
  num_individuals <- length(unique(qopt_data$Order))
  seq_len_n <- seq_len(num_individuals)
  
  # Replicate each sequence element 20 times (for 20 runs)
  seq_replicated <- rep(seq_len_n, each = 20)
  
  # Assign the replicated sequence as the SampleID
  qopt_data$SampleID <- factor(seq_replicated, levels = seq_len_n)  # Factor with ordered levels
  
  # Create the bar plot
  plot_data <- qopt_data %>%
    pivot_longer(cols = starts_with("Cluster"), names_to = "Cluster", values_to = "Proportion") %>%
    arrange(Order)  # Arrange the data by the "Order" column
  
  ggplot(plot_data, aes(x = SampleID, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Sample ID", y = "Admixture Proportion") +
    theme(axis.text.x = element_text(size = 6)) +
    ggtitle(paste("Admixture Plot for K =", K))
}

# Plot for a specific K value, e.g., K=4
plot_admixture(1)
plot_admixture(2)
plot_admixture(3)
plot_admixture(4)
plot_admixture(5)
