library(tidyverse)

### R code for pca 
setwd("/Users/natal/R/whales")

# read in the data
roh_data <- read.csv("plink_roh_prelim.csv")

# check the structure of the data
str(roh_data)

# Plink_ROH_association ---------------------------------------------------

# FROH calculation --------------------------------------------------------

# specify the total length of the genome in kilobases
total_genome_length_kb <- 3000000000  # replace with actual genome size

# calculate FROH for each individual (PLINK)
roh_data$plink_FROH <- roh_data$KB / total_genome_length_kb

# calculate FROH for each individual (PLINK)
roh_data$plink2_FROH <- roh_data$plink2_KB / total_genome_length_kb

# calculate FROH for each individual (bcf)
roh_data$bcf_FROH <- roh_data$bcf_KB / total_genome_length_kb

### filter for depth <5
# remove rows where depth < 5
filtered_data <- subset(roh_data, Depth >= 5)

# filter plink data
filtered_data_plink <- filtered_data[filtered_data$NSEG != 0 & filtered_data$KB != 0 & filtered_data$KBAVG != 0, ]

# filter plink2 data
filtered_data_plink2 <- filtered_data[filtered_data$plink2_NSEG != 0 & filtered_data$plink2_KB != 0 & filtered_data$plink2_KBAVG != 0, ]

# filter bcftools data
filtered_data_bcf <- filtered_data[filtered_data$bcf_NSEG != 0 & filtered_data$bcf_KB != 0 & filtered_data$bcf_KBAVG != 0, ]

#### PLINK1 associations
# linear regression: association between depth and nseg
nseg_model_plink1 <- lm(NSEG ~ Depth, data = filtered_data_plink)
summary(nseg_model_plink1)

# linear regression: association between depth and KB
kb_model_plink1 <- lm(KB ~ Depth, data = filtered_data_plink)
summary(kb_model_plink1)

# linear regression: association between depth and KBAVG
kbavg_model_plink1 <- lm(KBAVG ~ Depth, data = filtered_data_plink)
summary(kbavg_model_plink1)

# linear regression: association between depth and FROH
froh_model_plink1 <- lm(plink_FROH ~ Depth, data = filtered_data_plink)
summary(froh_model_plink1)

#### PLINK2 associations
# linear regression: association between depth and nseg
nseg_model_plink2 <- lm(plink2_NSEG ~ Depth, data = filtered_data_plink2)
summary(nseg_model_plink2)

# linear regression: association between depth and KB
kb_model_plink2 <- lm(plink2_KB ~ Depth, data = filtered_data_plink2)
summary(kb_model_plink2)

# linear regression: association between depth and KBAVG
kbavg_model_plink2 <- lm(plink2_KBAVG ~ Depth, data = filtered_data_plink2)
summary(kbavg_model_plink2)

# linear regression: association between depth and FROH
froh_model_plink2 <- lm(plink2_FROH ~ Depth, data = filtered_data_plink2)
summary(froh_model_plink2)

#### bcf associations
# linear regression: association between depth and nseg
nseg_model_bcf <- lm(bcf_NSEG ~ Depth, data = filtered_data_bcf)
summary(nseg_model_bcf)

# linear regression: association between depth and KB
kb_model_bcf <- lm(bcf_KB ~ Depth, data = filtered_data_bcf)
summary(kb_model_bcf)

# linear regression: association between depth and KBAVG
kbavg_model_bcf <- lm(bcf_KBAVG ~ Depth, data = filtered_data_bcf)
summary(kbavg_model_bcf)

# linear regression: association between depth and FROH
froh_model_bcf <- lm(bcf_FROH ~ Depth, data = filtered_data_bcf)
summary(froh_model_bcf)


# Load necessary library
library(tidyverse)

# Set working directory
setwd("/Users/natal/R/whales")

# Read in the data
roh_data <- read.csv("plink_roh_prelim.csv")

# Specify the total length of the genome in kilobases
total_genome_length_kb <- 3000000000  # replace with actual genome size in kilobases

# Calculate FROH for each individual
roh_data$plink_FROH <- roh_data$KB / total_genome_length_kb
roh_data$plink2_FROH <- roh_data$plink2_KB / total_genome_length_kb
roh_data$bcf_FROH <- roh_data$bcf_KB / total_genome_length_kb

# Filter for depth >= 5
filtered_data <- subset(roh_data, Depth >= 5)
filtered_data_plink <- filtered_data[filtered_data$NSEG != 0 & filtered_data$KB != 0 & filtered_data$KBAVG != 0, ]
filtered_data_plink2 <- filtered_data[filtered_data$plink2_NSEG != 0 & filtered_data$plink2_KB != 0 & filtered_data$plink2_KBAVG != 0, ]
filtered_data_bcf <- filtered_data[filtered_data$bcf_NSEG != 0 & filtered_data$bcf_KB != 0 & filtered_data$bcf_KBAVG != 0, ]

# Define models and labels
models <- list(
  plink_NSEG = lm(NSEG ~ Depth, data = filtered_data_plink),
  plink_KB = lm(KB ~ Depth, data = filtered_data_plink),
  plink_KBAVG = lm(KBAVG ~ Depth, data = filtered_data_plink),
  plink_FROH = lm(plink_FROH ~ Depth, data = filtered_data_plink),
  
  plink2_NSEG = lm(plink2_NSEG ~ Depth, data = filtered_data_plink2),
  plink2_KB = lm(plink2_KB ~ Depth, data = filtered_data_plink2),
  plink2_KBAVG = lm(plink2_KBAVG ~ Depth, data = filtered_data_plink2),
  plink2_FROH = lm(plink2_FROH ~ Depth, data = filtered_data_plink2),
  
  bcf_NSEG = lm(bcf_NSEG ~ Depth, data = filtered_data_bcf),
  bcf_KB = lm(bcf_KB ~ Depth, data = filtered_data_bcf),
  bcf_KBAVG = lm(bcf_KBAVG ~ Depth, data = filtered_data_bcf),
  bcf_FROH = lm(bcf_FROH ~ Depth, data = filtered_data_bcf)
)

# Extract p-values for the Depth coefficient from each model
p_values <- sapply(models, function(model) {
  summary(model)$coefficients["Depth", "Pr(>|t|)"]
})

# Create a data frame and add an asterisk next to p-values < 0.05, with formatting
p_values_table <- data.frame(
  Model = names(p_values),
  P_Value = ifelse(
    p_values < 0.05, 
    paste0(formatC(p_values, format = "f", digits = 4), "*"), 
    formatC(p_values, format = "f", digits = 4)
  )
)

# Print the p-values table
print(p_values_table)

