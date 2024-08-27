# Load the necessary libraries
library(ggplot2)
library(dplyr)

# Create a data frame with the heterozygosity values and species labels
data <- data.frame(
  species = factor(c("Blue whale", "Fin whale", "Humpback whale", 
                     "Minke whale", "N.At. right whale", "Sei whale A", "Sei whale B", 
                     "Gray whale", "Eastern gray whale", "Western gray whale"),
                   levels = c("Blue whale", "Fin whale", "Humpback whale", 
                              "Minke whale", "N.At. right whale", "Sei whale A", "Sei whale B", 
                              "Gray whale", "Eastern gray whale", "Western gray whale")),
  heterozygosity = c(0.00295, 0.00185, 0.0012, 
                     0.0006, 0.00175, 0.00075, 0.00077, 
                     0.0007, 0.000462, 0.000483)
)

# Create the bar plot with different colors for each species
ggplot(data, aes(x = species, y = heterozygosity, fill = species)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Species", y = "Heterozygosity", title = "Heterozygosity of Different Whale Species") +
  scale_fill_manual(values = c("Blue whale" = "skyblue", 
                               "Fin whale" = "#316d16", 
                               "Gray whale" = "#cc79a7", 
                               "Humpback whale" = "#d55e00", 
                               "Minke whale" = "#f0e442", 
                               "N.At. right whale" = "#9882cf", 
                               "Sei whale A" = "turquoise", 
                               "Sei whale B" = "#b87333", 
                               "Eastern gray whale" = "#ac0c06", 
                               "Western gray whale" = "#16537e"))


