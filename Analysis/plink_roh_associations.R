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

# calculate FROH for each individual (bcf)
roh_data$bcf_FROH <- roh_data$bcf_KB / total_genome_length_kb

#### more stringent PLINK parameters; no association above certain depth, but only finding large ROHs

### individuals with failed roh calculation are excluded
# remove rows where NSEG, KB, or KBAVG are 0
filtered_data <- roh_data[roh_data$NSEG != 0 & roh_data$KB != 0 & roh_data$KBAVG != 0, ]

# check the structure of the filtered data
str(filtered_data)

# linear regression: association between depth and nseg
nseg_model <- lm(NSEG ~ Depth, data = filtered_data)
summary(nseg_model)

# plot for nseg vs depth
plot(filtered_data$Depth, filtered_data$NSEG, 
     main = "Association between Depth and NSEG", 
     xlab = "Depth", ylab = "NSEG")
abline(nseg_model, col = "blue")
text(filtered_data$Depth, filtered_data$NSEG, labels = filtered_data$Ind, pos = 4, cex = 0.7, col = "gray40")

# linear regression: association between depth and KB
kb_model <- lm(KB ~ Depth, data = filtered_data)
summary(kb_model)

# plot for KB vs depth
plot(filtered_data$Depth, filtered_data$KB, 
     main = "Association between Depth and KB", 
     xlab = "Depth", ylab = "KB")
abline(kb_model, col = "blue")
text(filtered_data$Depth, filtered_data$KB, labels = filtered_data$Ind, pos = 4, cex = 0.7, col = "gray40")

# linear regression: association between depth and KBAVG
kbavg_model <- lm(KBAVG ~ Depth, data = filtered_data)
summary(kbavg_model)

# plot for KBAVG vs depth
plot(filtered_data$Depth, filtered_data$KBAVG, 
     main = "Association between Depth and KBAVG", 
     xlab = "Depth", ylab = "KBAVG")
abline(kbavg_model, col = "blue")
text(filtered_data$Depth, filtered_data$KBAVG, labels = filtered_data$Ind, pos = 4, cex = 0.7, col = "gray40")

### compare west vs east 

# count the number of entries for each POP
pop_counts <- table(filtered_data$POP)
print(pop_counts)

# calculate mean for each POP
mean_values <- aggregate(cbind(NSEG, KB, KBAVG) ~ POP, data = filtered_data, FUN = mean)
print(mean_values)

# t-test for NSEG
nseg_west <- filtered_data$NSEG[filtered_data$POP == "west"]
nseg_east <- filtered_data$NSEG[filtered_data$POP == "east"]
nseg_ttest <- t.test(nseg_west, nseg_east)
print(nseg_ttest)

# t-test for KB
kb_west <- filtered_data$KB[filtered_data$POP == "west"]
kb_east <- filtered_data$KB[filtered_data$POP == "east"]
kb_ttest <- t.test(kb_west, kb_east)
print(kb_ttest)

# t-test for KBAVG
kbavg_west <- filtered_data$KBAVG[filtered_data$POP == "west"]
kbavg_east <- filtered_data$KBAVG[filtered_data$POP == "east"]
kbavg_ttest <- t.test(kbavg_west, kbavg_east)
print(kbavg_ttest)


# plink2_association ------------------------------------------------------

### individuals with failed roh calculation are excluded
# remove rows where plink2_NSEG, plink2_KB, or plink2_KBAVG are 0
filtered_data <- roh_data[roh_data$plink2_NSEG != 0 & roh_data$plink2_KB != 0 & roh_data$plink2_KBAVG != 0, ]
# remove rows where depth < 5
filtered_data <- subset(filtered_data, Depth >= 7)

# check the structure of the filtered data
str(filtered_data)

# linear regression: association between depth and plink2_NSEG
plink2_NSEG_model <- lm(plink2_NSEG ~ Depth, data = filtered_data)
summary(plink2_NSEG_model)

# plot for plink2_NSEG vs depth
plot(filtered_data$Depth, filtered_data$plink2_NSEG, 
     main = "Association between Depth and plink2_NSEG", 
     xlab = "Depth", ylab = "plink2_NSEG")
abline(plink2_NSEG_model, col = "blue")
text(filtered_data$Depth, filtered_data$plink2_NSEG, labels = filtered_data$Ind, pos = 4, cex = 0.7, col = "gray40")

# linear regression: association between depth and plink2_KB
plink2_KB_model <- lm(plink2_KB ~ Depth, data = filtered_data)
summary(plink2_KB_model)

# plot for plink2_KB vs depth
plot(filtered_data$Depth, filtered_data$plink2_KB, 
     main = "Association between Depth and plink2_KB", 
     xlab = "Depth", ylab = "plink2_KB")
abline(plink2_KB_model, col = "blue")
text(filtered_data$Depth, filtered_data$plink2_KB, labels = filtered_data$Ind, pos = 4, cex = 0.7, col = "gray40")

# linear regression: association between depth and plink2_KBAVG
plink2_KBavg_model <- lm(plink2_KBAVG ~ Depth, data = filtered_data)
summary(plink2_KBavg_model)

# plot for plink2_KBAVG vs depth
plot(filtered_data$Depth, filtered_data$plink2_KBAVG, 
     main = "Association between Depth and plink2_KBAVG", 
     xlab = "Depth", ylab = "plink2_KBAVG")
abline(plink2_KBavg_model, col = "blue")
text(filtered_data$Depth, filtered_data$plink2_KBAVG, labels = filtered_data$Ind, pos = 4, cex = 0.7, col = "gray40")

### compare west vs east 

# count the number of entries for each POP
pop_counts <- table(filtered_data$POP)
print(pop_counts)

# calculate mean for each POP
mean_values <- aggregate(cbind(plink2_NSEG, plink2_KB, plink2_KBAVG) ~ POP, data = filtered_data, FUN = mean)
print(mean_values)

# t-test for plink2_NSEG
plink2_NSEG_west <- filtered_data$plink2_NSEG[filtered_data$POP == "west"]
plink2_NSEG_east <- filtered_data$plink2_NSEG[filtered_data$POP == "east"]
plink2_NSEG_ttest <- t.test(plink2_NSEG_west, plink2_NSEG_east)
print(plink2_NSEG_ttest)

# t-test for plink2_KB
plink2_KB_west <- filtered_data$plink2_KB[filtered_data$POP == "west"]
plink2_KB_east <- filtered_data$plink2_KB[filtered_data$POP == "east"]
plink2_KB_ttest <- t.test(plink2_KB_west, plink2_KB_east)
print(plink2_KB_ttest)

# t-test for plink2_KBAVG
plink2_KBavg_west <- filtered_data$plink2_KBAVG[filtered_data$POP == "west"]
plink2_KBavg_east <- filtered_data$plink2_KBAVG[filtered_data$POP == "east"]
plink2_KBavg_ttest <- t.test(plink2_KBavg_west, plink2_KBavg_east)
print(plink2_KBavg_ttest)


# Bcftools_ROH_association ------------------------------------------------

# remove rows where bcf_bcf_NSEG, bcf_KB, or bcf_KBAVG are 0
filtered_data <- roh_data[roh_data$bcf_NSEG != 0 & roh_data$bcf_KB != 0 & roh_data$bcf_KBAVG != 0, ]

# check the structure of the filtered data
str(filtered_data)

# linear regression: association between depth and bcf_NSEG
bcf_NSEG_model <- lm(bcf_NSEG ~ Depth, data = filtered_data)
summary(bcf_NSEG_model)

# plot for bcf_NSEG vs depth
plot(filtered_data$Depth, filtered_data$bcf_NSEG, 
     main = "Association between Depth and bcf_NSEG", 
     xlab = "Depth", ylab = "bcf_NSEG")
abline(bcf_NSEG_model, col = "blue")
text(filtered_data$Depth, filtered_data$bcf_NSEG, labels = filtered_data$Ind, pos = 4, cex = 0.7, col = "gray40")

# linear regression: association between depth and bcf_KB
bcf_KB_model <- lm(bcf_KB ~ Depth, data = filtered_data)
summary(bcf_KB_model)

# plot for bcf_KB vs depth
plot(filtered_data$Depth, filtered_data$bcf_KB, 
     main = "Association between Depth and bcf_KB", 
     xlab = "Depth", ylab = "bcf_KB")
abline(bcf_KB_model, col = "blue")
text(filtered_data$Depth, filtered_data$bcf_KB, labels = filtered_data$Ind, pos = 4, cex = 0.7, col = "gray40")

# linear regression: association between depth and bcf_KBAVG
bcf_KBavg_model <- lm(bcf_KBAVG ~ Depth, data = filtered_data)
summary(bcf_KBavg_model)

# plot for bcf_KBAVG vs depth
plot(filtered_data$Depth, filtered_data$bcf_KBAVG, 
     main = "Association between Depth and bcf_KBAVG", 
     xlab = "Depth", ylab = "bcf_KBAVG")
abline(bcf_KBavg_model, col = "blue")
text(filtered_data$Depth, filtered_data$bcf_KBAVG, labels = filtered_data$Ind, pos = 4, cex = 0.7, col = "gray40")

### compare west vs east 

# count the number of entries for each POP
pop_counts <- table(filtered_data$POP)
print(pop_counts)

# calculate mean for each POP
mean_values <- aggregate(cbind(bcf_NSEG, bcf_KB, bcf_KBAVG) ~ POP, data = filtered_data, FUN = mean)
print(mean_values)

# t-test for bcf_NSEG
bcf_NSEG_west <- filtered_data$bcf_NSEG[filtered_data$POP == "west"]
bcf_NSEG_east <- filtered_data$bcf_NSEG[filtered_data$POP == "east"]
bcf_NSEG_ttest <- t.test(bcf_NSEG_west, bcf_NSEG_east)
print(bcf_NSEG_ttest)

# t-test for bcf_KB
bcf_KB_west <- filtered_data$bcf_KB[filtered_data$POP == "west"]
bcf_KB_east <- filtered_data$bcf_KB[filtered_data$POP == "east"]
bcf_KB_ttest <- t.test(bcf_KB_west, bcf_KB_east)
print(bcf_KB_ttest)

# t-test for bcf_KBAVG
bcf_KBavg_west <- filtered_data$bcf_KBAVG[filtered_data$POP == "west"]
bcf_KBavg_east <- filtered_data$bcf_KBAVG[filtered_data$POP == "east"]
bcf_KBavg_ttest <- t.test(bcf_KBavg_west, bcf_KBavg_east)
print(bcf_KBavg_ttest)


# FROH calculation --------------------------------------------------------

# specify the total length of the genome in kilobases
total_genome_length_kb <- 3000000000  # replace with actual genome size

# calculate FROH for each individual (PLINK)
roh_data$plink_FROH <- roh_data$KB / total_genome_length_kb

# calculate FROH for each individual (bcf)
roh_data$bcf_FROH <- roh_data$bcf_KB / total_genome_length_kb

# inspect the results
head(roh_data[c("Ind", "FROH")])
View(roh_data)

write.csv(roh_data, "/Users/natal/Documents/Purdue/Whale Project/ROHs/ROH_data.csv")


# plot results ------------------------------------------------------------

library(tidyverse)

glimpse(roh_data)

# plink plot kb vs nseg 
ggplot(data = roh_data) +
  aes(x = KB, y = NSEG, color = POP) +
  geom_point()

# bcftools plot kb vs nseg 
ggplot(data = roh_data) +
  aes(x = bcf_KB, y = bcf_NSEG, color = POP) +
  geom_point()

# plink plot kb vs froh 
ggplot(data = roh_data) +
  aes(x = plink_FROH, y = NSEG, color = POP) +
  geom_point()

# bcftools plot kb vs froh
ggplot(data = roh_data) +
  aes(x = bcf_FROH, y = bcf_NSEG, color = POP) +
  geom_point()

# # plink plot kb vs froh
# ggplot(data = roh_data) +
#   aes(x = plink_FROH, color = POP) +
#   geom_point()
# 
# # bcftools plot kb vs froh 
# ggplot(data = roh_data) +
#   aes(x = bcf_FROH, color = POP) +
#   geom_boxplot() + 
#   geom_jitter()

filtered_data <- roh_data[roh_data$NSEG != 0 & roh_data$KB != 0 & roh_data$KBAVG != 0, ]

# plink froh boxplot
ggplot(data = filtered_data) +
  aes(x = POP, y = plink_FROH, color = POP) +
  geom_boxplot() +  # hide outliers to avoid overlap with jitter
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # jitter points for better visibility
  labs(x = "Population", y = "FROH") +
  theme_minimal()

# bcftools froh boxplot
ggplot(data = roh_data) +
  aes(x = POP, y = bcf_FROH, color = POP) +
  geom_boxplot() +  # hide outliers to avoid overlap with jitter
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # jitter points for better visibility
  labs(x = "Population", y = "FROH") +
  theme_minimal()









