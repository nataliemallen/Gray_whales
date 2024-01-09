setwd("/Users/natal/R")
kin<-read.csv("nobo_relate.csv",header = T)
#check number of sites compared and remove bad sites. In here sites less than 9950000 are removed
hist(kin$nSites)

#what is this for?? use or no?
#kin <- kin[-which(kin$nSites < 16800000),]
dev.off()
# population assignment

pop<-read.csv("/Users/natal/R/nobo_pop.csv", header = T)

head(pop)
#colnames(pop)[1]<-"ind"


#change the a and b values of kin data set to population labels
library(DataCombine)
r1 <- data.frame(from = pop$ind,
                 to=pop$pop)
r2 <- data.frame(from = pop$ind,
                 to=pop$sample_ID)
colnames(kin)[1]<-"a"
kin$a<-as.character(kin$a)
kin_new <- FindReplace(data = kin, Var = "a", replaceData = r1,
                       from = "from", to = "to", exact = T)

head(kin_new)
kin_new$b<-as.character(kin_new$b)

kin_new <- FindReplace(data = kin_new, Var = "b", replaceData = r1,
                       from = "from", to = "to", exact = T)

#change the ind_a and ind_b values of kin data set to sample IDs
kin_new$a<-as.character(kin_new$a)
kin_new <- FindReplace(data = kin_new, Var = "a", replaceData = r2,
                       from = "from", to = "to", exact = T)
kin_new$b<-as.character(kin_new$b)
kin_new <- FindReplace(data = kin_new, Var = "b", replaceData = r2,
                       from = "from", to = "to", exact = T)


east <- kin_new[which(kin_new$a =="east"  & kin_new$b =="east" ),]
west <- kin_new[which(kin_new$a =="west"  & kin_new$b =="west" ),]
unknown<- kin_new[which(kin_new$a =="unknown"  & kin_new$b =="unknown" ),]

#colors for plots
library(RColorBrewer)
transparent_blue <- rgb(0, 0, 1, alpha=0.5)

transparent_red <- rgb(1, 0, 0, alpha=0.5)

transparent_yel <- rgb(0, 1, 0, alpha=0.5)

transparent_pur <- rgb(1, 0, 1, alpha=0.5)

transparent_gre <- rgb(0, 1, 1, alpha=0.5)

#KING vs R1 for all pop 
plot(east$R1, east$KING, pch=20,cex=1, col="red", xlim = c(0.1,1),ylim=c(-0.2,0.5),
     xlab="R1", ylab="KING")+points(west$R1, west$KING, pch=20,cex=1, col="blue")+points(unknown$R1, unknown$KING, pch=20,cex=1, col="forestgreen")

polygon(c(0.423128,1,1,0.423128),c(0.1767767,0.1767767,0.3535534,0.3535534),col=transparent_color,
        border = NA) # Parent Offspring
#polygon(c(0.330533,0.427657,0.427657,0.330533),c(0.1767767,0.1767767,0.3535534,0.3535534),col=light_pal[4],
#border = NA) # Full Sib

polygon(c(0.2,0.5,0.5,0.2),c(0.08838835,0.08838835,0.1767767,0.1767767),col=transparent_red,
        border = NA) # Half Sib
polygon(c( 0.07142857,0.4210526,0.4210526,0.07142857),c(0.04419417,0.04419417,0.08838835,0.08838835),col=transparent_yel,
        border = NA) # First cousins
polygon(c(0,0.4,0.4,0),c(0.04419417,0.04419417,-0.04949,-0.04949),col=transparent_pur,
        border = NA) # Unrelated

legend(c(0.7,0.7),c(0.1,0), text.font = 1,
       legend=c("1st degree","2nd degree","3rd degree", "Unrelated"),
       fill = c(transparent_blue, transparent_red, transparent_yel, transparent_pur), bty='n', border = NA,
       pt.cex=2, x.intersp=1, text.width=2, xpd=T,y.intersp=1.5)

#KING vs R1 for pop east
plot(east$R1, east$KING, pch=20,cex=1, col="red", xlim = c(0.1,1),ylim=c(-0.2,0.5),
     xlab="R1", ylab="KING")

polygon(c(0.423128,1,1,0.423128),c(0.1767767,0.1767767,0.3535534,0.3535534),col=transparent_blue,
        border = NA) # Parent Offspring
polygon(c(0.330533,0.427657,0.427657,0.330533),c(0.1767767,0.1767767,0.3535534,0.3535534),col=transparent_gre,
border = NA) # Full Sib
polygon(c(0.2,0.5,0.5,0.2),c(0.08838835,0.08838835,0.1767767,0.1767767),col=transparent_red,
        border = NA) # Half Sib
polygon(c( 0.07142857,0.4210526,0.4210526,0.07142857),c(0.04419417,0.04419417,0.08838835,0.08838835),col=transparent_yel,
        border = NA) # First cousins
polygon(c(0,0.4,0.4,0),c(0.04419417,0.04419417,-0.04949,-0.04949),col=transparent_pur,
        border = NA) # Unrelated

legend(c(0.7,0.7),c(0.1,0), text.font = 1,
       legend=c("1st degree","2nd degree","3rd degree", "Unrelated"),
       fill = c(transparent_blue, transparent_red, transparent_yel, transparent_pur), bty='n', border = NA,
       pt.cex=2, x.intersp=1, text.width=2, xpd=T,y.intersp=1.5)

#KING vs R1 for pop west
plot(west$R1, west$KING, pch=20,cex=1, col="blue", xlim = c(0.1,1),ylim=c(-0.2,0.5),
     xlab="R1", ylab="KING")

polygon(c(0.423128,1,1,0.423128),c(0.1767767,0.1767767,0.3535534,0.3535534),col=transparent_blue,
        border = NA) # Parent Offspring
#polygon(c(0.330533,0.427657,0.427657,0.330533),c(0.1767767,0.1767767,0.3535534,0.3535534),col=light_pal[4],
#border = NA) # Full Sib
polygon(c(0.2,0.5,0.5,0.2),c(0.08838835,0.08838835,0.1767767,0.1767767),col=transparent_red,
        border = NA) # Half Sib
polygon(c( 0.07142857,0.4210526,0.4210526,0.07142857),c(0.04419417,0.04419417,0.08838835,0.08838835),col=transparent_yel,
        border = NA) # First cousins
polygon(c(0,0.4,0.4,0),c(0.04419417,0.04419417,-0.04949,-0.04949),col=transparent_pur,
        border = NA) # Unrelated

legend(c(0.7,0.7),c(0.1,0), text.font = 1,
       legend=c("1st degree","2nd degree","3rd degree", "Unrelated"),
       fill = c(transparent_blue, transparent_red, transparent_yel, transparent_pur), bty='n', border = NA,
       pt.cex=2, x.intersp=1, text.width=2, xpd=T,y.intersp=1.5)

#R0 vs R1 for all pop 
plot(east$R1, east$R0, pch=20,cex=1, col="red", xlim = c(0.1,1),ylim=c(-0.2,1.0),
     xlab="R1", ylab="KING")+points(west$R1, west$R0, pch=20,cex=1, col="blue")+points(unknown$R1, unknown$R0, pch=20,cex=1, col="forestgreen")

#R0 vs R1 for pop east
plot(east$R1, east$R0, pch=20,cex=1, col="red", xlim = c(0.1,1),ylim=c(-0.2,1.0),
     xlab="R1", ylab="R0")

#R0 vs R1 for pop west
plot(west$R1, west$R0, pch=20,cex=1, col="blue", xlim = c(0.1,1),ylim=c(-0.2,1.0),
     xlab="R1", ylab="R0")
