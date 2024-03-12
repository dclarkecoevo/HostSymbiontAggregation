#Code for seperating and cleaning data for host-parasite 
#  aggregation -- Fall 2019 -- Poecilia Reticulata-Gyrodactylus Spp

#Libraries needed for analysis and seperation 
library(plyr)
library(vegan)
library(tidyverse)
library(ggplot2)
library(lmodel2)
library(car)
library(lme4)
require(visreg)
require(arm)
require(lmodel2)
library(readr)


##Borrowed bits from Dr. Jessica Stephenson's code -- See ecology Guppy paper 


#Working Directory
getwd()#Show you current working directory
setwd("~/Desktop/DataforGuppy-GyroAgg")# Change directory to where the file is stored that you downloaded frokm github
#Import new dataset to work with in R
#### Load in data
dframe1 <- read_csv("Data/FINAL_TT1_New.csv", 
                      col_types = cols(Length = col_number(), 
                                       Wgt = col_number(), site2 = col_character(), 
                                       year = col_character()))

dframe1 <- read_csv("Data/FINAL_TT1_NEW12.csv", 
                    col_types = cols(Length = col_number(), 
                                     Wgt = col_number(), site2 = col_character(), 
                                     year = col_character()))


#eliminating any NA from dataset
dframe1<-filter(dframe1, Length != ".") #Eliminates the NA for length
dframe1<-filter(dframe1, Wgt != ".") #Eliminates the NA for wgt
dframe1<-filter(dframe1, drainage2 != ".")
dframe1<-filter(dframe1, course3 != ".")
dframe1<-filter(dframe1, river2 != ".")
dframe1<-filter(dframe1, Length > 0.0001)
dframe1<-filter(dframe1, Wgt > 0.0001)


###### Check the data is being read correctly
summary(dframe1)

##### Make sure variables are read as either numeric or factors - this also gives a description of each of the columns in the datafile
dframe1$Wgt <- as.numeric(dframe1$Wgt)# fish weight
dframe1$Length<-as.numeric(dframe1$Length)# fish Length
dframe1$gyro<-as.numeric(dframe1$gyro) # count of gyrodactylid parasites on that fish
dframe1$year<-as.factor(dframe1$year) # year the fish was sampled
dframe1$sex<-as.character(dframe1$sex)
dframe1$site2<-as.factor(dframe1$site2)
dframe1$drainage2<-as.factor(dframe1$drainage2)
dframe1$river<-as.factor(dframe1$river2)
dframe1$course3<-as.factor(dframe1$course3)
dframe1$presG<-as.factor(dframe1$presG)
dframe1$river2<-as.factor(dframe1$river2)
dframe1$season<-as.factor(dframe1$season)




dframe1$site2<-as.factor(dframe1$site2) # arbitrary site number within the course of each river - not useful alone but necessary for the next step:

###### Create nested spatial variables
dframe1 <- within(dframe1, {nriver<-factor(drainage2:river2)
ncourse <- factor(drainage2:river2:course3)
nsite <- factor(drainage2:river2:course3:site2)
nsiteyear<- factor(drainage2:river2:course3:site2:year)})



###############################NEW CODE################################
#Seperate males and females

dfm <- filter(dframe1, sex=="m")
dff<-filter(dframe1, sex=="f")
dfj<-filter(dframe1, sex=="j")


############ standardised scaled mass index - include the lake for this bit

lmodel2(log(as.numeric(Wgt)+1)~log(as.numeric(Length)+1), data=dframe1)
lmodel2(log(as.numeric(Wgt)+1)~log(as.numeric(Length)+1), data=dff)
lmodel2(log(as.numeric(Wgt)+1)~log(as.numeric(Length)+1), data=dfm)
lmodel2(log(as.numeric(Wgt)+1)~log(as.numeric(Length)+1), data=dfj)
### use the SMA slope

summary(dframe1$Length) ### use the mean as L0
summary(dff$Length)
summary(dfm$Length)
summary(dfj$Length)

dframe1$SMI<-(dframe1$Wgt)*((15.2/dframe1$Length)^0.2644078)
dff$SMI<-(dff$Wgt)*((15.00/dff$Length)^0.3727303)
dfm$SMI<-(dfm$Wgt)*((15.74/dfm$Length)^0.3136675)
dfj$SMI<-(dfj$Wgt)*((10.18/dfj$Length)^0.09381226)

###### scale the SMI15 variable so that males and females both have a mean of 0
dff<-subset(dframe1, sex=="f")

################### AND MALES
dfm<-subset(dframe1, sex=="m")
dfm$scSMI<-(dfm$SMI-mean(dfm$SMI))

###### scale the JUV

dfj$scSMI<-(dfj$SMI-mean(dfj$SMI))

############## and so that course within sex also have a mean of 0
dfflow<-subset(dff, course3=="lower")

dffup<-subset(dff, course3=="upper")

############## and so that course within sex also have a mean of 0
dfmlow<-subset(dfm, course3=="lower")

dfmup<-subset(dfm, course3=="upper")


#### Calculating for Upper and Lower
#Upper

lmodel2(log(as.numeric(Wgt)+1)~log(as.numeric(Length)+1), data=dffup)
summary(dffup$Length)
dffup$SMI<-(dffup$Wgt)*((18.32/dffup$Length)^0.4350118)

lmodel2(log(as.numeric(Wgt)+1)~log(as.numeric(Length)+1), data=dfmup)
summary(dfmup$Length)
dfmup$SMI<-(dfmup$Wgt)*((16.31/dfmup$Length)^0.2972246)


#Lower

lmodel2(log(as.numeric(Wgt)+1)~log(as.numeric(Length)+1), data=dfflow)
summary(dfflow$Length)
dfflow$SMI<-(dfflow$Wgt)*((17.00/dfflow$Length)^0.336538)


lmodel2(log(as.numeric(Wgt)+1)~log(as.numeric(Length)+1), data=dfmlow)
summary(dfmlow$Length)
dfmlow$SMI<-(dfmlow$Wgt)*((15.0/dfmlow$Length)^0.3307505)


#Calculating the mean and standard deviation standardized
dfmup$scSMI2<-(dfmup$SMI-mean(dfmup$SMI))
dfmup$normSMI<-(dfmup$scSMI2-mean(dfmup$scSMI2))

dfmlow$scSMI2<-(dfmlow$SMI-mean(dfmlow$SMI))
dfmlow$normSMI<-(dfmlow$scSMI2-mean(dfmlow$scSMI2))

dffup$scSMI2<-(dffup$SMI-mean(dffup$SMI))
dffup$normSMI<-(dffup$scSMI2-mean(dffup$scSMI2))


dfflow$scSMI2<-(dfflow$SMI-mean(dfflow$SMI))
dfflow$normSMI<-(dfflow$scSMI2-mean(dfflow$scSMI2))



#Creating the new dataframes
dfmc<-rbind(dfmlow, dfmup)
dffc<-rbind(dfflow,dffup)



dframe2<-rbind(dffc,dfmc)

#Getting the H and P for all 
Hall<-ddply(dframe2,.(nsiteyear), nrow)   #Takes the total rows in the dataframe  to calculate total H 

Pall<-ddply(dframe2,.(nsiteyear), summarize,Ptot=sum(gyro), MeanP=mean(gyro),VarP=var(gyro),meanLen=mean(Length),varLen=var(Length), meanW=mean(Wgt), varW=var(Wgt), sdlen=sd(Length), minlen=min(Length), maxlen=max(Length), meanBC=mean(scSMI2), varBC=var(scSMI2))  #Takes the sum of gyro vector in the data frame to calculate total P, mean P, and variance in P 

HPdataall<-cbind(Hall,Ptot=Pall$Ptot,MeanP=Pall$MeanP,VarP=Pall$VarP,meanLen=Pall$meanLen,varLen=Pall$varLen, meanW=Pall$meanW, varW=Pall$varW, sdlen=Pall$sdlen,minlen=Pall$minlen,maxlen=Pall$maxlen, meanBC=Pall$meanBC, varBC=Pall$varBC)  #Binds the two estimates together with the site/river/sample name to give you the information together 
##Getting Total H and P for females

Hf<-ddply(dffc,.(nsiteyear), nrow)   #Takes the total rows in the dataframe  to calculate total H 

Pf<-ddply(dffc,.(nsiteyear), summarize, Ptot=sum(gyro), MeanP=mean(gyro),VarP=var(gyro),meanLen=mean(Length),varLen=var(Length), meanW=mean(Wgt), varW=var(Wgt), sdlen=sd(Length), minlen=min(Length), maxlen=max(Length), meanBC=mean(scSMI2), varBC=var(scSMI2))  #Takes the sum of gyro vector in the data frame to calculate total P, mean P, and variance in P 

HPdataf<-cbind(Hf,Ptot=Pf$Ptot,MeanP=Pf$MeanP,VarP=Pf$VarP,meanLen=Pf$meanLen,varLen=Pf$varLen, meanW=Pf$meanW, varW=Pf$varW, sdlen=Pf$sdlen, minlen=Pf$minlen, maxlen=Pf$maxlen, meanBC =Pf$meanBC, varBC=Pf$varBC)  #Binds the two estimates together with the site/river/sample name to give you the information together 
##Getting total H and P for males

Hm<-ddply(dfmc,.(nsiteyear), nrow)   #Takes the total rows in the dataframe  to calculate total H 

Pm<-ddply(dfmc,.(nsiteyear),summarize,Ptot=sum(gyro), MeanP=mean(gyro),VarP=var(gyro),meanLen=mean(Length),varLen=var(Length), meanW=mean(Wgt), varW=var(Wgt),sdLen=sd(Length),minlen=min(Length), maxlen=max(Length), meanBC=mean(scSMI2), varBC=var(scSMI2))  #Takes the sum of gyro vector in the data frame to calculate total P, mean P, and variance in P 

HPdatam<-cbind(Hm,Ptot=Pm$Ptot,MeanP=Pm$MeanP,VarP=Pm$VarP,meanLen=Pm$meanLen,varLen=Pm$varLen, meanW=Pm$meanW, varW=Pm$varW, sdlen = Pm$sdLen,minlen=Pm$minlen,maxlen=Pm$maxlen, meanBC=Pm$meanBC, varBC=Pm$varBC)  #Binds the two estimates together with the site/river/sample name to give you the information together 



#Getting total H and P for jouviniles and females

Hfj<-ddply(dfjc,.(nsiteyear), nrow)   #Takes the total rows in the dataframe  to calculate total H 

Pfj<-ddply(dfjc,.(nsiteyear),summarize,Ptot=sum(gyro), MeanP=mean(gyro),VarP=var(gyro),meanLen=mean(Length),varLen=var(Length), meanW=mean(Wgt), varW=var(Wgt),sdLen=sd(Length),minlen=min(Length), maxlen=max(Length), meanBC=mean(scSMI2), varBC=var(scSMI2))  #Takes the sum of gyro vector in the data frame to calculate total P, mean P, and variance in P 

HPdatafj<-cbind(Hfj,Ptot=Pfj$Ptot,MeanP=Pfj$MeanP,VarP=Pfj$VarP,meanLen=Pfj$meanLen,varLen=Pfj$varLen, meanW=Pfj$meanW, varW=Pfj$varW, sdlen = Pfj$sdLen,maxlen=Pfj$maxlen, SMI=Pfj$scSMI2, varBC=Pf$varBC)  #Binds the two estimates together with the site/river/sample name to give you the information together 

HPdatafj<-drop_na(HPdatafj)

#Convert to data frame -- Females
HPdff<-data.frame(Site=HPdataf$nsite,
                  H=HPdataf$V1,
                  P=HPdataf$Ptot, 
                  MeanP=HPdataf$MeanP,
                  VarP=HPdataf$VarP, 
                  MeanLen=HPdataf$meanLen, 
                  VarLen = HPdataf$varLen, 
                  MeanW=HPdataf$meanW, 
                  VarW=HPdataf$varW, 
                  sdlen=HPdataf$sdlen,
                  minlen=HPdataf$minlen,
                  maxlen=HPdataf$maxlen,
                  meanBC=HPdataf$meanBC,
                  varBC=HPdataf$varBC
                  
)


#Convert to data frame -- Males
HPdfm<-data.frame(Site=HPdatam$nsite,
                  H=HPdatam$V1,
                  P=HPdatam$Ptot, 
                  MeanP=HPdatam$MeanP,
                  VarP=HPdatam$VarP, 
                  MeanLen=HPdatam$meanLen, 
                  VarLen = HPdatam$varLen, 
                  MeanW=HPdatam$meanW, 
                  VarW=HPdatam$varW,
                  sdlen=HPdatam$sdlen,
                  minlen=HPdatam$minlen,
                  maxlen=HPdatam$maxlen,
                  meanBC=HPdatam$meanBC,
                  varBC=HPdatam$varBC
                  
)


#Convery to data frame -- all
HPdfall<-data.frame(Site=HPdataall$nsite,
                    H=HPdataall$V1,P=HPdataall$Ptot,
                    MeanP=HPdataall$MeanP,VarP=HPdataall$VarP, 
                    MeanLen=HPdataall$meanLen, 
                    VarLen = HPdataall$varLen, 
                    MeanW=HPdataall$meanW, 
                    VarW=HPdataall$varW,
                    sdlen=HPdataall$sdlen, 
                    minlen=HPdataall$minlen,
                    maxlen=HPdataall$maxlen,
                    meanBC=HPdataall$meanBC,
                    varBC=HPdataall$varBC
)

#Convery to data frame -- fj
HPdffj<-data.frame(Site=HPdatafj$nsite,
                   H=HPdatafj$V1,
                   P=HPdatafj$Ptot,
                   MeanP=HPdatafj$MeanP,
                   VarP=HPdatafj$VarP, 
                   MeanLen=HPdatafj$meanLen, 
                   VarLen = HPdatafj$varLen, 
                   MeanW=HPdatafj$meanW, 
                   VarW=HPdatafj$varW,
                   sdlen=HPdatafj$sdlen, 
                   prev=HPdatafj$prev,
                   minlen=HPdatafj$minlen,
                   maxlen=HPdatafj$maxlen
)
#----------------
# Now that we have total estimates we subset out any populations who did not have gyrodactylus 
#----------------
#Females
HPdff= filter(HPdff, P > 4) #taking out populations with less than 10 gyrodactylus (H and P limitations to be determined later)
HPdff=filter(HPdff, H > 3) # taking out populations with less than 10 hosts


#Males
HPdfm= filter(HPdfm, P > 4) #taking out populations with less than 10 gyrodactylus (H and P limitations to be determined later)
HPdfm=filter(HPdfm, H > 3) # taking out populations with less than 10 hosts


#all
HPdfall= filter(HPdfall, P > 4) #taking out populations with less than 10 gyrodactylus (H and P limitations to be determined later)
HPdfall=filter(HPdfall, H > 3) # taking out populations with less than 10 hosts

#fj
HPdffj= filter(HPdffj, P > 4) #taking out populations with less than 10 gyrodactylus (H and P limitations to be determined later)
HPdffj=filter(HPdffj, H > 3) # taking out populations with less than 10 hosts

######Host het






#Save this data frame for Python Analysis 
write.csv(HPdff, "HP_Aggregation_Trinidad_females100.csv")
write.csv(HPdfm, "HP_Aggregation_Trinidad_males100.csv")
write.csv(HPdfall, "HP_Aggregation_Trinidad_all100.csv")
write.csv(HPdffj, "HP_Aggregation_Trinidad_Female&Juv1.csv")




###########

#Extra

#data frames to see the distribtion of parasites among size classes

quantile(dframe_f$Length)
dframe_s<-filter(dframe_f, length<92)
dframe_m<-filter(dframe_f, length>92 & length<202)
dframe_l<-filter(dframe_f, length>202)
hist(dframe_s$gyro, breaks=100,xlim = range(0,100), ylim=range(0,250))
hist(dframe_m$gyro,breaks=100,xlim = range(0,100), ylim=range(0,250))
hist(dframe_l$gyro, breaks=100,xlim = range(0,100), ylim=range(0,250))

dframe_f_u<-filter(dframe_f, course=="upper")
quantile(dframe_f_u$length)
dframe_f_l<-filter(dframe_f, course=="lower")
quantile(dframe_f_l$length)


smfu<-filter(dframe_f, length<=92)
smfu<-filter(smf, course=="upper")
mfu<-filter(dframe_f,length>=92 & length<=202)
mfu<-filter(mf, course=="upper")
lfu<-filter(dframe_f,length<=202)
lfu<-filter(lf, course=="upper")


smfl<-filter(smf, course=="lower")

mfl<-filter(mf, course=="lower")

lfl<-filter(lf, course=="lower")
par(mfrow=c(3,1))

hist(smfu$gyro, breaks=100 )
hist(mfu$gyro, breaks=100)
hist(lfu$gyro, breaks=100)
hist(smfl$gyro)
hist(mfl$gyro)
hist(lfl$gyro)

dframej<-filter(dframe1, class=="j")
hist(dframej$gyro)

#attaining total H and P from each river/site/year combination
plot(minlen~jonvar, FS_trinf)
plot(maxlen~jonvar, FS_trinf)



###Calculating a sd(length)~Mean(length) heterogeneity
lenhetf<-resid(lm(sdlen~MeanLen, data=HPdff))
lenhetm<-resid(lm(sdlen~MeanLen, data=HPdfm))
lenhetall<-resid(lm(sdlen~MeanLen, data=HPdfall))
lenhetfj<-resid(lm(sdlen~MeanLen, data=HPdffj))
HPdfm$sdlenhet<-lenhetm
HPdff$sdlenhet<-lenhetf
HPdfall$sdlenhet<-lenhetall
HPdffj$sdlenhet<-lenhetfj
