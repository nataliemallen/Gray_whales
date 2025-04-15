### R code for pca 
setwd("/Users/natal/R/whales")

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(ggsignif)

# Read the CSV file
data <- read.csv("heterozygosity.csv")
#View(data)

# Reorder ID based on heterozygosity
data$ID <- factor(data$ID, levels = data$ID[order(data$heterozygosity, decreasing = TRUE)])

str(data)

# Plot
ggplot(data, aes(x = ID, y = heterozygosity, fill = POP)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "heterozygosity by ID",
       x = "ID",
       y = "heterozygosity") +
  scale_fill_manual(values = c("east" = "#1E88E5", "west" = "#D81B60")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# new plot 
# Ensure the 'POP' column is treated as a factor with correct order
data$POP <- factor(data$POP, levels = c("east", "west"))

# Define the colors for each POP
POP_colors <- c("east" = "#1E88E5",  
                 "west" = "#e51e8a")

filtered_data <- data %>%
  group_by(POP) %>%
  filter(n() > 2)

p <- ggplot(filtered_data, aes(x = POP, y = heterozygosity, color = POP)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(x = POP, y = heterozygosity), width = 0.2, alpha = 0.6) +
  theme_minimal() +
  ylim(0, 0.0016) +
  scale_color_manual(values = POP_colors) +
  
  # Customize x-axis and y-axis labels
  labs(x = "Sample Site", y = "Heterozygosity", color = "Sample Site") +  # Change labels
  
  # Increase the size of x-axis labels
  theme(axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12))

print(p)

# perform the wilcoxon rank-sum test
wilcox_test <- wilcox.test(heterozygosity ~ POP, data = data)

# print the results
print(wilcox_test)



