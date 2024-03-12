#This script is to run alongside the AggregationSIMay2021_V4.Rmd file 
#to print the graphs we want to use for publication
cbp1 <- c( "#E69F00", "#56B4E9")
cbp2<-c("#F0E442")
cbp3<-c("#CC79A7")
library(ggplot2)
library(gridExtra)

#Limpet panel code
#Limpet boxplot of diffvar
Limpgraph1<-ggplot(Limp, aes(y=diffvar, x=L, fill=L))+
geom_boxplot(aes(),  width=0.5, outlier.shape=NA)+
theme_classic()+
geom_hline(yintercept=0, linetype="dashed")+
xlab("Limpet Populations")+
  ylab(" ")+
#ylab("Difference in variance (Obs.-Exp.)")+
scale_fill_manual(values=cbp2)+
theme_classic()+
geom_jitter(color="grey33")+
  theme(legend.position="none")+
  annotate("text",size=10,label="a", fontface= 1)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=21), axis.ticks.x=element_blank(),axis.text.x=element_blank())
#Limpet Partial resid diffvar versus mean length
Limpgraph2<-ggplot(VSrsL2, aes(MeanLRS, visregRes))+
  geom_smooth(method="lm")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_jitter(data=VSrsL2, aes(MeanLRS, visregRes))+
  xlab("Scaled Mean Length")+
  ylab("")+
  theme_classic()+
  theme(legend.position="none")+
  annotate("text", size=10,label="b", fontface= 1)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=21))
#Limpet Partial resid diffvar versus variance in length
Limpgraph3<-ggplot(VSrsL, aes(VarLRS, visregRes))+
  geom_smooth(method="lm")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_jitter(data=VSrsL, aes(VarLRS, visregRes))+
  xlab("Scaled Variance in Length")+
  ylab(" ")+
  theme_classic()+
  theme(legend.position="none")+
  annotate("text",  size=10,label="c", fontface= 1)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=21))
  
#Panel of limpet figures 
m<-ggplot()+geom_blank()

grid.arrange(m,Limpgraph1,Limpgraph2,Limpgraph3, ncol=2, nrow=2)
ggsave("Limpetgraphs.pdf")

#pdf("Limpetgraphs_final.pdf")
#Limpetgraphs<-grid.arrange(m,Limpgraph1,Limpgraph2,Limpgraph3, ncol=2, nrow=2)
#print(Limpetgraphs)
#dev.off()
#Guppy aggregation panel
#Guppy boxplot for diffvar
GuppyAgg1<-ggplot(FS_trin, aes(G,diffvar, fill=G))+
  geom_boxplot( width=0.5,outlier.shape=NA)+
  theme_classic()+
  geom_hline(yintercept=0, linetype="dashed")+
  xlab("Guppy Populations")+
  ylab(" ")+
  #ylab("Difference in variance (Obs.-Exp.)")+
  scale_fill_manual(values=cbp3)+
  geom_jitter(width=0.4)+
  theme(legend.position="none")+
  ylim(-0.6,0.5)+
  annotate("text", y=.5, x=.5, size=10,label="a", fontface= 1)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=22))
  
#Guppy Partial resid diffvar plot versus Mean Length

GuppyAgg2<-ggplot(VSDLr, aes(MeanLenRS, visregRes))+
  geom_smooth(method="lm")+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_jitter(data=VSDLr, aes(MeanLenRS, visregRes))+
  xlab("Scaled Mean Length (mm)")+
  ylab("")+
  theme_classic()+
  ylim(-0.6,0.5)+
  annotate("text", y=.5, x=-3, size=10,label="b", fontface= 1)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=20))
  

#Guppy Partial resid diffvar plot versus Sex
GuppyAgg3<-ggplot(VSrs, aes(sex, visregRes, fill=sex))+
  geom_boxplot()+
  geom_hline(yintercept=0,linetype="dashed")+
  geom_jitter(data=VSrs, aes(sex, visregRes))+
  xlab("Sex")+
  ylab("")+
  theme_classic()+
  scale_fill_manual(values=cbp1)+
  theme(legend.position="none")+
  ylim(-0.6,0.5)+
  annotate("text", y=.5, x=.5, size=10,label="c", fontface= 1)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=20))
  

#Panel of guppy figures
grid.arrange(m,GuppyAgg1,GuppyAgg2,GuppyAgg3, nrow=2, ncol=2)

#If we wanted to see all six diffvar plots together 
#grid.arrange(Limpgraph1,Limpgraph2,Limpgraph3,GuppyAgg1,GuppyAgg3,GuppyAgg2, nrow=2)

#Prevalence and Max parasite ratio graphs

#Prevalence graphs
GuppyPrev1<-ggplot(FS_trin, aes(G,DiffPrev, fill=G))+
  geom_boxplot(width=0.5,outlier.shape=NA)+
  theme_classic()+
  geom_hline(yintercept=0, linetype="dashed")+
  xlab(" ")+
  ylab("Difference in Prevalence (Obs.-Exp.)")+
  scale_fill_manual(values=cbp3)+
  geom_jitter(width=0.4)+
  theme(legend.position="none")+
  ylim(-0.6,0.65)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=18),plot.margin=margin(10,10,10,10))+
  annotate("text", y=.6, x=.50, size=10,label="a", fontface= 1)

#-0.2,0.5
#diffprev partial resid vs mean length
GuppyPrev2<-ggplot(VSPpr, aes(MeanLenRS, visregRes))+
  theme_classic()+
  geom_smooth(method="lm")+
  geom_jitter(data=VSPpr, aes(MeanLenRS, visregRes))+
  xlab(" ")+
  ylab(" ")+
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(-0.6,0.65)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=21),plot.margin=margin(10,10,10,10))+
  annotate("text", y=.6, x=-2.7, size=10,label="b", fontface= 1)

#diffprev partial resid vs sex
GuppyPrev3<-ggplot(VSPpr2, aes(sex, visregRes,))+
  theme_classic()+
  geom_boxplot(aes(fill=sex),outlier.shape=NA)+
  geom_jitter(data=VSPpr2, aes(sex, visregRes))+
  xlab(" ")+
  ylab(" ")+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_fill_manual(values=cbp1)+
  ylim(-0.6,0.65)+
  annotate("text", y=.6, x=.5, size=10, label="c", fontface= 1,)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=21), legend.position="none",plot.margin=margin(10,10,10,10) )
  


#Prevalence only panel
#grid.arrange(GuppyPrev1,GuppyPrev2,GuppyPrev3, nrow=1)

#Max parasite ratio graphs
#Logged max parasite ratio 
GuppyMPR1<-ggplot(FS_trin, aes(G,log10(DiffPRatio),fill=G))+
  geom_boxplot( width=0.5,outlier.shape=NA)+
  theme_classic()+
  geom_hline(yintercept=0, linetype="dashed")+
  xlab("Guppy Populations")+
  ylab("Log Ratio of Max parasite  (Obs./ Exp.)")+
  scale_fill_manual(values=cbp3)+
  geom_jitter(width=0.4)+
  theme(legend.position="none")+
  ylim(-0.6,0.65)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=21),plot.margin=margin(10,10,10,10))+
  annotate("text", y=.6, x=.5, size=10,label="d", fontface= 1)

GuppyMPR2<-ggplot(VRMLr, aes(MeanLenRS, visregRes))+
  geom_smooth(method="lm")+
  geom_jitter(data=VRMLr, aes(MeanLenRS,visregRes))+
  xlab("Scaled Mean Length (mm)")+
  ylab(" ")+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_color_manual(values=cbp1)+
  theme_classic()+
  ylim(-0.6,0.65)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=21),plot.margin=margin(10,10,10,10) )+
  annotate("text", y=.6, x=-3, size=10,label="e", fontface= 1)
  
  
GuppyMPR3<-ggplot(VSMPr, aes(sex, visregRes, fill=sex))+
  geom_boxplot()+
  geom_jitter(data=VSMPr, aes(sex,visregRes))+
  xlab("Sex")+
  ylab(" ")+
  geom_hline(yintercept=0, linetype="dashed")+
  scale_fill_manual(values=cbp1)+
  theme_classic()+
  theme(legend.position="none")+
  ylim(-0.6,0.65)+
  theme(axis.text=element_text(size=25),axis.title=element_text(size=21),plot.margin=margin(10,10,10,10))+
  annotate("text", y=.6, x=.5, size=10,label="f", fontface= 1)


#Max Parasite ratio panel only
#grid.arrange(GuppyMPR1,GuppyMPR2,GuppyMPR3, nrow=1)

grid.arrange(GuppyPrev1,GuppyPrev2,GuppyPrev3, GuppyMPR1,GuppyMPR2,GuppyMPR3, nrow=2,  heights=c(2,2))
