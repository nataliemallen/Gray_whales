### R code for admixture plot
setwd("/Users/natal/R/whales")

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Function to process the data for a specific K value
process_data <- function(K) {
  # Load .csv dataframe
  file_name <- paste0("whale_k", K, "_combined.csv")
  qopt_data <- read.csv(file_name, header = TRUE)
  
  # Convert cluster columns to numeric
  cluster_columns <- paste0("Cluster", 1:K)
  qopt_data[cluster_columns] <- lapply(qopt_data[cluster_columns], as.numeric)
  
  # Replace "east" with "E" and "west" with "W"
  qopt_data$POP <- ifelse(qopt_data$POP == "east", "E", "W")
  
  # Order the rows by the "Order" column
  qopt_data <- qopt_data[order(qopt_data$Order), ]
  
  # Calculate the average proportion for each cluster for each individual
  averaged_data <- qopt_data %>%
    group_by(Order, POP) %>%
    summarise(across(all_of(cluster_columns), mean)) %>%
    ungroup()
  
  # Add a column for K
  averaged_data$K <- K
  
  return(averaged_data)
}

# Combine data for all K values
combined_data <- bind_rows(lapply(2:5, process_data))

# Reshape the data to long format for plotting
plot_data <- combined_data %>%
  pivot_longer(cols = starts_with("Cluster"), 
               names_to = "Cluster", 
               values_to = "Proportion") %>%
  arrange(K, Order)

# Create the bar plot
ggplot(plot_data, aes(x = factor(Order), y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ K, scales = "free_x") +
  labs(x = "Individual", y = "Admixture Proportion") +
  theme(axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank()) +
  ggtitle("Admixture Plots for K=2, 3, 4, and 5")
