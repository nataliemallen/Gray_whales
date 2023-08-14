###from Andrew Black

library(readxl)
library(ggplot2)
whale_metadata <- read_excel("Library/CloudStorage/Box-Box/Personal/Postdoc_Purdue/whale_popgen/whale_metadata.xlsx")

#Breadth, post filtering
ggplot(data=whale_metadata, aes(y=breadth_post, x=reorder(ID,breadth_post),fill=pop))+geom_bar(stat="identity")+theme_classic()+scale_fill_manual("", values =c("west"="blue","east"="black"))+xlab("Sample (N=73)")+ylab("Genome Breadth")+theme( axis.ticks.x=element_blank(),axis.text.x = element_blank())+ geom_hline(yintercept=50, linetype="dashed", color = "white", linewidth=1)+ylim(0,100)+theme(legend.position="top")
#Depth, post filtering
ggplot(data=whale_metadata, aes(y=depth_post, x=reorder(ID,depth_post),fill=pop))+geom_bar(stat="identity")+theme_classic()+scale_fill_manual("", values =c("west"="blue","east"="black"))+xlab("Sample (N=73)")+ylab("Mean Depth Of Coverage")+theme( axis.ticks.x=element_blank(),axis.text.x = element_blank())+ geom_hline(yintercept=2, linetype="dashed", color = "white", linewidth=1)+ylim(0,12)+theme(legend.position="top")
