##Graphics for P. reticulata aggregation 

#This file contains scirpt for generating extra graphs for Guppy-Gyrodactylus script. 

#Packages
library(plyr)
library(vegan)
library(tidyverse)
library(ggplot2)
library(lmodel2)
library(car)
library(lme4)
require(visreg)
require(arm)


dframe1 <- read.csv("Data/Trini_data_wCond.csv")    # Select "Trini_data_wCond.csv"

#eliminating any NA from dataset
dframe1<-filter(dframe1, length != ".") #Eliminates the NA for length
dframe1<-filter(dframe1, wgt != ".") #Eliminates the NA for wgt
#resetting pres/absence for prevelance

dframe1$presG[dframe1$presG=="1"]<-"0"
dframe1$presG[dframe1$presG=="2"]<-"1"

###### Check the data is being read correctly
summary(dframe1)

##### Make sure variables are read as either numeric or factors - this also gives a description of each of the columns in the datafile
dframe1$W <- as.numeric(dframe1$wgt) # fish weight
dframe1$gyro<-as.numeric(dframe1$gyro) # count of gyrodactylid parasites on that fish
dframe1$year<-as.factor(dframe1$year) # year the fish was sampled
dframe1$Cond<-as.numeric(dframe1$Cond)
dframe1$length<-as.numeric(dframe1$length)

dframe1$presG<-as.numeric(dframe1$presG) # presence/absence of gyrodactylid parasites
dframe1$prestrich<-as.factor(dframe1$prestrich) # presence/absence of trichodina
dframe1$presWS<-as.factor(dframe1$presWS) # presence/absence white spot
dframe1$presApio<-as.factor(dframe1$presApio) # presence/absence apiosoma
dframe1$presfung<-as.factor(dframe1$presfung) # presence/absence fungal infection 
dframe1$presDig<-as.factor(dframe1$presDig) # presence/absence of digenean metacercariae 
dframe1$presCN<-as.factor(dframe1$presCN) # presence/absence of camallanus nematodes 

dframe1$richness.no.gyro<-as.numeric(dframe1$richness.no.gyro) # richness score not including gyrodactylid parasites
dframe1$richness.w.gyro<-as.numeric(dframe1$richness.w.gyro) # richness score including gyrodactylid parasites

dframe1$site<-as.factor(dframe1$site) # arbitrary site number within the course of each river - not useful alone but necessary for the next step:

###### Create nested spatial variables
dframe1 <- within(dframe1, {nriver<-factor(drainage:river)
    ncourse <- factor(drainage:river:course)
    nsite <- factor(drainage:river:course:site)
    nsiteyear<- factor(drainage:river:course:site:year)})

dframe1f<-filter(dframe1, class== "f")
dframe1j<-filter(dframe1, class=="j")
dframe1fj<-rbind(dframe1f,dframe1j)
quantile(dframe1f$length)

dframe1s<-filter(dframe1f, length < 92)
dframe1m<-filter(dframe1f, length > 92 & length < 202)
dframe1l<-filter(dframe1f, length > 202)
#Looking at a SSAD of different size classes

#Juveniles 
Agg_Hist_Juve=ggplot(data=dframe1j, aes(x=gyro, y=..density..))+
  geom_histogram( aes(x=gyro, y=..density..), fill="dark green")+
  theme_classic()+
  ylab("Density")+
  xlab("Number of Gyrodactylus")+
  theme(text = element_text(size=35))+
  ylim(0,1.1)
  
pdf("Figures/Agg_Hist_Juve.pdf") #Saving the figure as a pdf
Agg_Hist_Juve
dev.off()

#Small 
Agg_Hist_Small=ggplot(dframe1s, aes(gyro))+
    geom_histogram(aes(x=gyro, y=..density..), fill="dark blue")+
  theme_classic()+
  ylab("Density")+
  xlab("Number of Gyrodactylus")+
  theme(text = element_text(size=35))+
  ylim(0,1)
  
pdf("Figures/Agg_Hist_Small.pdf") #Saving the figure as a pdf
Agg_Hist_Small
dev.off()

ggplot(dframe1m, aes(gyro))+
  geom_histogram(aes(x=gyro, y=..density..), fill="dark red")+
  theme_classic()+ylab("Density")+
  xlab("Number of Gyrodactylus")+
  theme(text = element_text(size=35))+
  ylim(0,1)
  
pdf("Figures/Agg_TPL_all.pdf") #Saving the figure as a pdf
Agg_TPL_all
dev.off()

ggplot(dframe1l, aes(gyro))+
  geom_histogram(aes(x=gyro, y=..density..), fill="dark gray")+
  theme_classic()+
  ylab("Density")+
  xlab("Number of Gyrodactylus")+
  theme(text = element_text(size=35))+
  ylim(0,1)
  #scale_x_continuous(limits=c(0,60))
  
summary(dframe1$gyro)

#Histogram for looking at every population 

hist(dframe1$gyro[nsiteyear=="yarra:yarra:upper:1:2006"], main="yarra:yarra:upper:1:2006", breaks=20)
text(2,30,"diff= 0.076") 


hist(dframe1$gyro[nsiteyear=="cunupia:dyke:lower:1:2006"], main="cunupia:dyke:lower:1:2006", breaks=20)
text(.5,10,"diff= -0.407") 


hist(dframe1$gyro[nsiteyear=="caroni:guanapo:lower:1:2004" ], main="caroni:guanapo:lower:1:2004", breaks=20)
text(.5,30,"diff= -0.363") 


hist(dframe1$gyro[nsiteyear=="caroni:aripo:upper:3:2006"], main="caroni:aripo:upper:3:2006", breaks=20)
text(1,20,"diff= -0.351") 


hist(dframe1$gyro[nsiteyear=="laseiva:laseiva:lower:1:2006"], main="laseiva:laseiva:lower:1:2006", breaks=20)
text(1,30,"diff= -0.252") 


hist(dframe1$gyro[nsiteyear=="caroni:lopinot:upper:1:2004"], main="caroni:lopinot:upper:1:2004", breaks=20)
text(1.5,15,"diff= -0.215")

hist(dframe1$gyro[nsiteyear=="marianne:marianne:lower:1:2004"], main="marianne:marianne:lower:1:2004", breaks=20)
text(3,1.5,"diff= 0.056") 
hist(dframe1$gyro[nsiteyear=="caroni:guanapo:lower:1:2003"], main="caroni:guanapo:lower:1:2003", breaks=20)
text(2,30,"diff= -0.142") 
hist(dframe1$gyro[nsiteyear=="guapo:guapo:lower:2:2006"], main="guapo:guapo:lower:2:2006", breaks=20)
text(3,2.5,"diff= -0.124") 
hist(dframe1$gyro[nsiteyear=="caroni:dyke:lower:3:2004"], main="caroni:dyke:lower:3:2004", breaks=20)
text(3,3,"diff= -0.100") 
hist(dframe1$gyro[nsiteyear=="marianne:marianne:upper:1:2004"], main="marianne:marianne:upper:1:2004", breaks=20)
text(1.5,10,"diff= -0.341")
hist(dframe1$gyro[nsiteyear=="caroni:caura:upper:1:2003"], main="caroni:caura:upper:1:2003", breaks=20)
text(1.5,15,"diff= -0.222")
hist(dframe1$gyro[nsiteyear=="yarra:yarra:upper:1:2004"], main="yarra:yarra:upper:1:2004", breaks=20)
text(3,6,"diff= -0.106")
hist(dframe1$gyro[nsiteyear=="caroni:guanapo:upper:1:2003"], main="caroni:guanapo:upper:1:2003", breaks=20)
text(1.5,20,"diff= -0.286")
hist(dframe1$gyro[nsiteyear=="marianne:marianne:upper:1:2003"], main="marianne:marianne:upper:1:2003", breaks=20)
text(2,40,"diff= -0.382")
hist(dframe1$gyro[nsiteyear=="oropuche:turure:lower:1:2004"], main="oropuche:turure:lower:1:2004", breaks=20)
text(2,30,"diff= -0.463")
hist(dframe1$gyro[nsiteyear=="caroni:caura:upper:3:2003"], main="caroni:caura:upper:3:2003", breaks=20)
text(2,30,"diff= -0.480")
hist(dframe1$gyro[nsiteyear=="caroni:arima:upper:1:2003"], main="caroni:arima:upper:1:2003", breaks=20)
text(2,30,"diff= -0.226")
hist(dframe1$gyro[nsiteyear=="matura:matura:lower:1:2006"], main="matura:matura:lower:1:2006", breaks=20)
text(2,20,"diff= -0.582")
hist(dframe1$gyro[nsiteyear=="caroni:arima:lower:2:2006"], main="caroni:arima:lower:2:2006", breaks=20)
text(20,15,"diff= -0.204")
hist(dframe1$gyro[nsiteyear=="caroni:aripo:upper:7:2008"], main="caroni:aripo:upper:7:2008", breaks=20)
text(3,10,"diff= -0.463")
hist(dframe1$gyro[nsiteyear=="caroni:aripo:upper:6:2003"], main="caroni:aripo:upper:6:2003", breaks=20)
text(2,20,"diff= -0.630")
hist(dframe1$gyro[nsiteyear=="oropuche:oropuche:lower:1:2006"], main="oropuche:oropuche:lower:1:2006", breaks=20)
text(2.5,10,"diff= -0.665")
hist(dframe1$gyro[nsiteyear=="caroni:caura:upper:1:2004"], main="caroni:caura:upper:1:2004", breaks=20)
text(3,30,"diff= -0.492")
hist(dframe1$gyro[nsiteyear=="caroni:dyke:lower:3:2003"], main="caroni:dyke:lower:3:2003", breaks=20)
text(2,30,"diff= -0.554")
hist(dframe1$gyro[nsiteyear=="yarra:yarra:lower:1:2003"], main="yarra:yarra:lower:1:2003", breaks=20)
text(6,20,"diff= -0.230")
hist(dframe1$gyro[nsiteyear=="caroni:aripo:lower:1:2003"], main="caroni:aripo:lower:1:2003", breaks=20)
text(2,30,"diff= -0.364")
hist(dframe1$gyro[nsiteyear=="caroni:aripo:upper:4:2006"], main="caroni:aripo:upper:4:2006", breaks=20)
text(4,20,"diff= -0.439")
hist(dframe1$gyro[nsiteyear=="caroni:arima:lower:2:2006"], main="caroni:arima:lower:2:2006", breaks=20)
text(20,10,"diff= -0.204")
hist(dframe1$gyro[nsiteyear=="caroni:aripo:upper:2:2008"], main="caroni:aripo:upper:2:2008", breaks=20)
text(10,10,"diff= -0.078")
hist(dframe1$gyro[nsiteyear=="caroni:lopinot:upper:1:2006"], main="caroni:lopinot:upper:1:2006", breaks=20)
text(7,15,"diff= -0.198")
hist(dframe1$gyro[nsiteyear=="caroni:lopinot:upper:1:2003"], main="caroni:lopinot:upper:1:2003", breaks=20)
text(2,30,"diff= -0.455")
hist(dframe1$gyro[nsiteyear=="caroni:guanapo:lower:1:2006"], main="caroni:guanapo:lower:1:2006", breaks=20)
text(3,10,"diff= -0.727")
hist(dframe1$gyro[nsiteyear=="marianne:marianne:lower:1:2003"], main="marianne:marianne:lower:1:2003", breaks=20)
text(7,30,"diff= -0.409")
hist(dframe1$gyro[nsiteyear=="caroni:dyke:lower:2:2006"], main="caroni:dyke:lower:2:2006", breaks=20)
text(30,30,"diff= 0.334")
hist(dframe1$gyro[nsiteyear=="caroni:aripo:lower:1:2008"], main="caroni:aripo:lower:1:2008", breaks=20)
text(20,10,"diff= 0.178")
hist(dframe1$gyro[nsiteyear=="mayaro:trib:lower:1:2006"], main="mayaro:trib:lower:1:2006", breaks=20)
text(40,6,"diff= 0.123")
hist(dframe1$gyro[nsiteyear=="caroni:aripo:upper:8:2003"], main="caroni:aripo:upper:8:2003", breaks=20)
text(5,30,"diff= -0.609")
hist(dframe1$gyro[nsiteyear=="caroni:dyke:lower:1:2006"], main="caroni:dyke:lower:1:2006", breaks=20)
text(30,2.0,"diff= -0.064")
hist(dframe1$gyro[nsiteyear=="caroni:arima:lower:2:2006"], main="caroni:arima:lower:2:2006", breaks=20)
text(15,15,"diff= -0.204")
hist(dframe1$gyro[nsiteyear=="caroni:lopinot:lower:1:2004"], main="caroni:lopinot:lower:1:2004", breaks=20)
text(20,30,"diff= -0.192")
hist(dframe1$gyro[nsiteyear=="caroni:lopinot:lower:1:2006"], main="caroni:lopinot:lower:1:2006", breaks=20)
text(10,15,"diff= -0.192")
hist(dframe1$gyro[nsiteyear=="caroni:dyke:lower:2:2003"], main="caroni:dyke:lower:2:2003", breaks=20)
text(15,30,"diff= -0.273")
hist(dframe1$gyro[nsiteyear=="caroni:aripo:lower:1:2006"], main="caroni:aripo:lower:1:2006", breaks=20)
text(20,10,"diff= -0.185")
hist(dframe1$gyro[nsiteyear=="caroni:maracas:upper:1:2003"], main="caroni:maracas:upper:1:2003", breaks=20)
text(20,30,"diff= -0.202")
hist(dframe1$gyro[nsiteyear=="caroni:lopinot:lower:1:2003"], main="caroni:lopinot:lower:1:2003", breaks=20)
text(50,50,"diff= 0.319")
#Extra subsetting for other plots

levels(dframe1$nsiteyear)

loop.vector <- 1:84
nsiteyear<-dframe1$nsiteyear
for (i in nsiteyear){
  pdf(file=a)
 a<-hist(dframe1$gyro[nsiteyear==i], main=paste(i))

}


plot_list = list()
for (i in nsiteyear) {
  hist(dframe1$gyro[nsiteyear==i], main=paste(i))
}

pdf("plots.pdf")

for (i in 1:84) {
  print(plot_list[[i]])
}
dev.off()




