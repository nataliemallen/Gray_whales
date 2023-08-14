### from Andrew Black

library(ggplot2)
meta <- read_excel("meta.xlsx")
View(meta)                  
cov<-as.matrix(read.table("~/prelim.cov"))
head(axes$values/sum(axes$values)*100)
[1] 5.153134 3.307909 2.903597
[4] 2.699114 2.535402 2.406617
PC1_3<-as.data.frame(axes$vectors[,1:3])
x<-cbind(PC1_3,meta)

ggplot(data=x, aes(y=V2, x=V1))+geom_point(size=7,pch=21,aes(fill=pop))+ theme_classic() + xlab("PC1 (5.15%)") +ylab("PC2 (3.31%)")+geom_hline(yintercept=0,linetype="dashed")+geom_vline(xintercept =0,linetype="dashed")+scale_fill_manual("Population", values=c("blue","black"))
