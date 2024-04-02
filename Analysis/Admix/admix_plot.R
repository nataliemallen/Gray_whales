### R code for k = 2 admixture plot
#.csv file columns: Sample_ID	Cluster1	Cluster2	Order	POP
#Clusters from .qopt; Order can be used to set order of bars 
setwd("/Users/natal/R")

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Load .csv dataframe
qopt_data <- read.csv("whale_k2_ordered.csv", header = TRUE)

# Sort the dataframe based on the "Order" column
qopt_data <- qopt_data[order(qopt_data$Order), ]

# Convert columns to numeric if needed
qopt_data$Cluster1 <- as.numeric(qopt_data$Cluster1)
qopt_data$Cluster2 <- as.numeric(qopt_data$Cluster2)

# Reshape the data to long format
plot_data <- qopt_data %>%
  pivot_longer(cols = c(Cluster1, Cluster2), names_to = "Cluster", values_to = "Proportion") %>%
  arrange(Order)  # Arrange the data by the "Order" column

# Create a sequence from 1 to n (number of rows in qopt_data)
seq_len_n <- seq_len(nrow(qopt_data))

# Replicate each sequence element twice (for two clusters)
seq_replicated <- rep(seq_len_n, each = 2)

# Assign the replicated sequence as the SampleID
plot_data$SampleID <- factor(seq_replicated, levels = seq_len_n)  # Factor with ordered levels

# Create the bar plot
ggplot(plot_data, aes(x = SampleID, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample ID", y = "Admixture Proportion")

ggplot(plot_data, aes(x = SampleID, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample ID", y = "Admixture Proportion") +
  theme(axis.text.x = element_text(size = 6))

### k = 3
# Load .csv dataframe
qopt_data <- read.csv("whale_k3_ordered.csv", header = TRUE)

# Sort the dataframe based on the "Order" column
qopt_data <- qopt_data[order(qopt_data$Order), ]

# Convert columns to numeric if needed
qopt_data$Cluster1 <- as.numeric(qopt_data$Cluster1)
qopt_data$Cluster2 <- as.numeric(qopt_data$Cluster2)
qopt_data$Cluster3 <- as.numeric(qopt_data$Cluster3)

# Reshape the data to long format
plot_data <- qopt_data %>%
  pivot_longer(cols = c(Cluster1, Cluster2, Cluster3), names_to = "Cluster", values_to = "Proportion") %>%
  arrange(Order)  # Arrange the data by the "Order" column

# Create a sequence from 1 to n (number of rows in qopt_data)
seq_len_n <- seq_len(nrow(qopt_data))

# Replicate each sequence element twice (for two clusters)
seq_replicated <- rep(seq_len_n, each = 3)

# Assign the replicated sequence as the SampleID
plot_data$SampleID <- factor(seq_replicated, levels = seq_len_n)  # Factor with ordered levels

# Create the bar plot
ggplot(plot_data, aes(x = SampleID, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample ID", y = "Admixture Proportion") +
  theme(axis.text.x = element_text(size = 6))


### k = 4
### R code for pca 
setwd("/Users/natal/R")

# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Load .csv dataframe
qopt_data <- read.csv("whale_k4_ordered.csv", header = TRUE)

# Sort the dataframe based on the "Order" column
qopt_data <- qopt_data[order(qopt_data$Order), ]

# Convert columns to numeric if needed
qopt_data$Cluster1 <- as.numeric(qopt_data$Cluster1)
qopt_data$Cluster2 <- as.numeric(qopt_data$Cluster2)
qopt_data$Cluster3 <- as.numeric(qopt_data$Cluster3)
qopt_data$Cluster4 <- as.numeric(qopt_data$Cluster4)

# Reshape the data to long format
plot_data <- qopt_data %>%
  pivot_longer(cols = c(Cluster1, Cluster2, Cluster3, Cluster4), names_to = "Cluster", values_to = "Proportion") %>%
  arrange(Order)  # Arrange the data by the "Order" column

# Create a sequence from 1 to n (number of rows in qopt_data)
seq_len_n <- seq_len(nrow(qopt_data))

# Replicate each sequence element twice (for two clusters)
seq_replicated <- rep(seq_len_n, each = 4)

# Assign the replicated sequence as the SampleID
plot_data$SampleID <- factor(seq_replicated, levels = seq_len_n)  # Factor with ordered levels

# Create the bar plot
ggplot(plot_data, aes(x = SampleID, y = Proportion, fill = Cluster)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample ID", y = "Admixture Proportion") +
  theme(axis.text.x = element_text(size = 6))
