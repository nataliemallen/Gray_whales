#Clear RStudio’s memory
rm(list = ls())

setwd("/Users/natal/R/whales")
kin<-read.csv("whale_nobo_relate.csv",header = T)
#check number of sites compared and remove bad sites. In here sites less than 9950000 are removed
hist(kin$nSites)

#what is this for?? use or no?
#kin <- kin[-which(kin$nSites < 16800000),]
dev.off()
# population assignment

pop<-read.csv("/Users/natal/R/whales/whale_nobo_pop.csv", header = T)

head(pop)
#colnames(pop)[1]<-"ind"

#change the a and b values of kin data set to population labels
library(DataCombine)
r1 <- data.frame(from = pop$ind,
                 to=pop$pop)
r2 <- data.frame(from = pop$ind,
                 to=pop$sample_ID)
colnames(kin)[1]<-"a"

kin_new <- kin
kin_new$ind_a <- kin_new$a #new
kin_new$ind_b <- kin_new$b #new

#replace with pop data
kin_new$a<-as.character(kin_new$a)
kin_new <- FindReplace(data = kin_new, Var = "a", replaceData = r1,
                       from = "from", to = "to", exact = T)
kin_new$b<-as.character(kin_new$b)
kin_new <- FindReplace(data = kin_new, Var = "b", replaceData = r1,
                       from = "from", to = "to", exact = T)

#change to individual names 
#View(kin_new)

#change the ind_a and ind_b values of kin data set to sample IDs
kin_new$ind_a<-as.character(kin_new$ind_a)
kin_new <- FindReplace(data = kin_new, Var = "ind_a", replaceData = r2,
                       from = "from", to = "to", exact = T)
kin_new$ind_b<-as.character(kin_new$ind_b)
kin_new <- FindReplace(data = kin_new, Var = "ind_b", replaceData = r2,
                       from = "from", to = "to", exact = T)

# Export kin_new
write.csv(kin_new, "kin_new.csv", row.names = FALSE)

all_individuals <- kin_new
east <- kin_new[which(kin_new$a =="east"  & kin_new$b =="east" ),]
west <- kin_new[which(kin_new$a =="west"  & kin_new$b =="west" ),]

###plot with dots as relationships - all individuals
# Set up the plot
plot(all_individuals$R1, all_individuals$KING, pch=20, cex=1, col="black", xlim=c(0.1,1), ylim=c(-0.2,0.5),
     xlab="R1", ylab="KING")


#bounds <- matrix(c(0.423128, Inf, 0.1767767, 0.3535534,   # Parent Offspring
                   #0.2, 0.5, 0.08838835, 0.1767767,       # Half Sib
                   #0.07142857, 0.4210526, 0.04419417, 0.08838835,  # First cousins
                   #-Inf, Inf, -Inf, 0.04419417),          # Unrelated
                # ncol = 4, byrow = TRUE)

#changed 1st degree bound to include parents and sibs
bounds <- matrix(c(0.38, Inf, 0.1767767, 0.3535534,   # Parent Offspring
                   0.2, 0.5, 0.08838835, 0.1767767,       # Half Sib
                   0.07142857, 0.4210526, 0.04419417, 0.08838835,  # First cousins
                   -Inf, Inf, -Inf, 0.04419417),          # Unrelated
                 ncol = 4, byrow = TRUE)

# Define colors based on the bounds
colors <- c("#E76BF3", "#00A5FF", "#00BC59", "#F8766D")

# Loop through each relationship category and plot the points with the corresponding color
for (i in 1:length(colors)) {
  points(kin_new$R1[kin_new$R1 >= bounds[i,1] & kin_new$R1 <= bounds[i,2] &
                      kin_new$KING >= bounds[i,3] & kin_new$KING <= bounds[i,4]], 
         kin_new$KING[kin_new$R1 >= bounds[i,1] & kin_new$R1 <= bounds[i,2] &
                        kin_new$KING >= bounds[i,3] & kin_new$KING <= bounds[i,4]], 
         pch=20, cex=1, col=colors[i])
}

# Add legend
legend(c(0.7, 0.7), c(0.1, 0), text.font=1,
       legend=c("1st degree", "2nd degree", "3rd degree", "Unrelated"),
       fill=colors, bty='n', border=NA,
       pt.cex=2, x.intersp=1, text.width=2, xpd=TRUE, y.intersp=1.5)


###plot with dots as relationships - east
# Set up the plot
plot(east$R1, east$KING, pch=20, cex=1, col="black", xlim=c(0.1,1), ylim=c(-0.2,0.5),
     xlab="R1", ylab="KING")


#bounds <- matrix(c(0.423128, Inf, 0.1767767, 0.3535534,   # Parent Offspring
                  # 0.2, 0.5, 0.08838835, 0.1767767,       # Half Sib
                  # 0.07142857, 0.4210526, 0.04419417, 0.08838835,  # First cousins
                  # -Inf, Inf, -Inf, 0.04419417),          # Unrelated
               #  ncol = 4, byrow = TRUE)

#changed 1st degree bound to include parents and sibs
bounds <- matrix(c(0.39, Inf, 0.1767767, 0.3535534,   # Parent Offspring
                   0.2, 0.5, 0.08838835, 0.1767767,       # Half Sib
                   0.07142857, 0.4210526, 0.04419417, 0.08838835,  # First cousins
                   -Inf, Inf, -Inf, 0.04419417),          # Unrelated
                 ncol = 4, byrow = TRUE)

# Define colors based on the bounds
colors <- c("#E76BF3", "#00A5FF", "#00BC59", "#F8766D")

# Loop through each relationship category and plot the points with the corresponding color
for (i in 1:length(colors)) {
  points(east$R1[east$R1 >= bounds[i,1] & east$R1 <= bounds[i,2] &
                     east$KING >= bounds[i,3] & east$KING <= bounds[i,4]], 
         east$KING[east$R1 >= bounds[i,1] & east$R1 <= bounds[i,2] &
                       east$KING >= bounds[i,3] & east$KING <= bounds[i,4]], 
         pch=20, cex=1, col=colors[i])
}

# Add legend
legend(c(0.7, 0.7), c(0.1, 0), text.font=1,
       legend=c("1st degree", "2nd degree", "3rd degree", "Unrelated"),
       fill=colors, bty='n', border=NA,
       pt.cex=2, x.intersp=1, text.width=2, xpd=TRUE, y.intersp=1.5)


###plot with dots as relationships - west
# Set up the plot
plot(west$R1, west$KING, pch=20, cex=1, col="black", xlim=c(0.1,1), ylim=c(-0.2,0.5),
     xlab="R1", ylab="KING")


#bounds <- matrix(c(0.423128, Inf, 0.1767767, 0.3535534,   # Parent Offspring
                  # 0.2, 0.5, 0.08838835, 0.1767767,       # Half Sib
                  # 0.07142857, 0.4210526, 0.04419417, 0.08838835,  # First cousins
                  # -Inf, Inf, -Inf, 0.04419417),          # Unrelated
                # ncol = 4, byrow = TRUE)

#changed 1st degree bound to include parents and sibs
bounds <- matrix(c(0.38, Inf, 0.1767767, 0.3535534,   # Parent Offspring
                   0.2, 0.5, 0.08838835, 0.1767767,       # Half Sib
                   0.07142857, 0.4210526, 0.04419417, 0.08838835,  # First cousins
                   -Inf, Inf, -Inf, 0.04419417),          # Unrelated
                 ncol = 4, byrow = TRUE)

# Define colors based on the bounds
colors <- c("#E76BF3", "#00A5FF", "#00BC59", "#F8766D")

# Loop through each relationship category and plot the points with the corresponding color
for (i in 1:length(colors)) {
  points(west$R1[west$R1 >= bounds[i,1] & west$R1 <= bounds[i,2] &
                     west$KING >= bounds[i,3] & west$KING <= bounds[i,4]], 
         west$KING[west$R1 >= bounds[i,1] & west$R1 <= bounds[i,2] &
                       west$KING >= bounds[i,3] & west$KING <= bounds[i,4]], 
         pch=20, cex=1, col=colors[i])
}

# Add legend
legend(c(0.7, 0.7), c(0.1, 0), text.font=1,
       legend=c("1st degree", "2nd degree", "3rd degree", "Unrelated"),
       fill=colors, bty='n', border=NA,
       pt.cex=2, x.intersp=1, text.width=2, xpd=TRUE, y.intersp=1.5)
