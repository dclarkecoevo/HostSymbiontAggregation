#Analysis for aggregation code. This step occurs after the python analysis where we calculate the feasible set. 

#This script is used for running linear model analysis of aggregation data from the Guppy-Gyrodactylus project. 
#Please make sure to change working directory to wherever you downloaded these files to ensure ease of code. 

rm(list = ls())

#Packages
#####
#Installing packages
#install.packages("vegan")
#install.packages("ggplot2")
#install.packages("lmodel2")
#install.packages("car")
#install.packages("lme4")
#install.packages("visreg")
#install.packages("arm")
#install.packages("tidyverse")
#install.packages("plyr")


#Loading Packages
library(vegan)
library(ggplot2)
library(lmodel2)
library(car)
library(lme4)
require(visreg)
require(arm)
library(tidyverse)
library(plyr)
library(corrplot)
source("HighstatLibV6.R")
library(DHARMa)
library(nlme)
#####


#Working Directory
getwd()#Show you current working directory
setwd("~/Desktop/Projects/DataforGuppy-GyroAgg")# Change directory to where the file is stored that you downloaded frokm github
#Import new dataset to work with in R
FS_trinf<-read.csv("Data/hostpara_tpl_data_females100.csv") #Select hostparasite_tpl_data_Trinidad_female
FS_trinm<-read.csv("Data/hostpara_tpl_data_males100.csv") #Select hostparasite_tpl_data_Trinidad_male


#Setting sex 
FS_trinf$sex<-"F"
FS_trinm$sex<-"M"




#Setting the correct class of variables
FS_trinf$fs_var_mean_log<-as.numeric(FS_trinf$fs_var_mean_log)
FS_trinm$fs_var_mean_log<-as.numeric(FS_trinm$fs_var_mean_log)


#Filter out courses we dont want
FS_trinf<-filter(FS_trinf, Course != "lake" ) 
FS_trinf<-filter(FS_trinf, Course != "mid" )#Filter out populations with hosts less than parasites
FS_trinm<-filter(FS_trinm, Course != "lake" ) 
FS_trinm<-filter(FS_trinm, Course != "mid" ) #Filter out populations with hosts less than parasites
FS_trinall<-filter(FS_trinall, Course != "lake" ) 
FS_trinall<-filter(FS_trinall, Course != "mid" )

#FS_trinf<-filter(FS_trinf, river != "Dyke")
#FS_trinall<-filter(FS_trinall, river != "Dyke")
#FS_trinm<-filter(FS_trinm, river != "Dyke")

FS_trinf<-filter(FS_trinf, river != "Guapo")
FS_trinall<-filter(FS_trinall, river != "Guapo")
FS_trinm<-filter(FS_trinm, river != "Guapo")

FS_trinf<-filter(FS_trinf, river != "Matura")
FS_trinall<-filter(FS_trinall, river != "Matura")
FS_trinm<-filter(FS_trinm, river != "Matura")

FS_trinf<-filter(FS_trinf, drainage != "Mayaro")
FS_trinall<-filter(FS_trinall, drainage != "Mayaro")
FS_trinm<-filter(FS_trinm, drainage != "Mayaro")

#variable setting

FS_trinf$year<-as.factor(FS_trinf$year) #Setting year as categorical
FS_trinm$year<-as.factor(FS_trinm$year) #Setting year as categorical
FS_trinall$year<-as.factor(FS_trinall$year) #Setting year as categorical


#Variance through Johnson method for length

jonvarf<-lm(VarLen~MeanLen, data=FS_trinf)
varresidf<-resid(jonvarf)
FS_trinf$jonvar<-varresidf

jonvarm<-lm(VarLen~MeanLen, data=FS_trinm)
varresidm<-resid(jonvarm)
FS_trinm$jonvar<-varresidm

jonvarall<-lm(VarLen~MeanLen, data=FS_trinall)
varresidall<-resid(jonvarall)
FS_trinall$jonvar<-varresidall

#Variance for Body condition
ResBCf<-lm(meanBC~MeanLen, data=FS_trinf)
FS_trinf$ResBC<-resid(ResBCf)
jonresf<-lm(varBC~ResBC, data=FS_trinf)
FS_trinf$jonres<-resid(jonresf)


ResBCm<-lm(meanBC~MeanLen, data=FS_trinm)
FS_trinm$ResBC<-resid(ResBCm)
jonresm<-lm(varBC~ResBC, data=FS_trinm)
FS_trinm$jonres<-resid(jonresm)

ResBCall<-lm(meanBC~MeanLen, data=FS_trinall)
FS_trinall$ResBC<-resid(ResBCall)
jonresall<-lm(varBC~ResBC, data=FS_trinall)
FS_trinall$jonres<-resid(jonresall)




#Calculating the unexplained variance by subtracting log FS variance from log observed variance
unexvarf=log10(FS_trinf$VarP)-log10(FS_trinf$fs_var) #Getting the difference between observed and FS variance
unexvarm=log10(FS_trinm$VarP)-log10(FS_trinm$fs_var) #Getting the difference between observed and FS variance
unexarall=log10(FS_trinall$VarP)-log10(FS_trinall$fs_var) #Getting the difference between observed and FS variance


FS_trinf$diffvar=unexvarf #Adding the difference to the dataset
FS_trinm$diffvar=unexvarm#Adding the difference to the dataset
FS_trinall$diffvar=unexarall#Adding the difference to the dataset

#Body Condition residuals 
FS_trinf$ResBC<-resid(lm(FS_trinf$meanBC~FS_trinf$MeanLen))
FS_trinf$ResBC<-as.numeric(scale(FS_trinf$ResBC))
FS_trinm$ResBC<-resid(lm(FS_trinm$meanBC~FS_trinm$MeanLen))
FS_trinm$ResBC<-as.numeric(scale(FS_trinm$ResBC))
FS_trinall$ResBC<-resid(lm(FS_trinall$meanBC~FS_trinall$MeanLen))
FS_trinall$ResBC<-as.numeric(scale(FS_trinall$ResBC))

#Calculating logMean
FS_trinf$logMean<-log10(FS_trinf$MeanP)
FS_trinm$logMean<-log10(FS_trinm$MeanP)
FS_trinall$logMean<-log10(FS_trinall$MeanP)



FS_trin<-rbind(FS_trinf,FS_trinm)

#Rescaling for model
#FS_trinall$jonres<-as.numeric(scale(FS_trinall$jonres))
#FS_trinf$jonres<-as.numeric(scale(FS_trinf$jonres))
FS_trin$jonresRS<-as.numeric(scale(FS_trin$jonres))
#FS_trinf$MeanLen<-as.numeric(scale(FS_trinf$MeanLen))
FS_trin$MeanLenRS<-as.numeric(scale(FS_trin$MeanLen))
FS_trin$MeanLenRS<-as.numeric(scale(FS_trin$MeanLen))
#Models for predicting which biological and ecological variables predict difference in variance
write.csv(FS_trin, "FS_trin.csv")
#Linear models for all with sex
lmallwsex<-lmer(diffvar~Course+ResBC+jonres+jonvar+MeanLen+logMean+sex+sex:MeanLen+sex:ResBC+(1|drainage),data=FS_trin)##Linear model for predicting  unexplained variance for all with multiple variables
summary(lmallwsex)## Summary for linear model 
vif(lmallwsex)#Checking  model for inflation factors
Anova(lmallwsex)#Looking at important of each predictor 
out<-capture.output(summary(lmallwsex))
cat("Results", out, file="Results/Summarylmallwsex.txt", sep="\n", append=FALSE)
out<-capture.output(Anova(lmallwsex))
cat("Results", out, file="Results/AnovaSummarylmallwsex.txt", sep="\n", append=FALSE)

#Some model validation
#~Course+ResBC+jonres+jonvar+MeanLen+logMean

#plot(lmallwsex)
#resid<-residuals(lmallwsex)
#hist(resid)
#plot(resid~FS_trin$year)
#plot(resid~FS_trin$Course)
#plot(resid~FS_trin$ResBC)
#plot(resid~FS_trin$jonres)
#plot(resid~FS_trin$jonvar)
#plot(resid~FS_trin$MeanLen)
#plot(resid~FS_trin$logMean)

#Linear models for females
lmf<-lm(diffvar~Course+ResBC+MeanLen+year+logMean+jonvar, data=FS_trinf) #Linear models for females
vif(lmf)#Checking female model for inflation factor
Anova(lmf) #Looking at important of each predictor 
summary(lmf) #Summary of linear model for females
out<-capture.output(summary(lmf))
cat("Results", out, file="Results/SummarylmFemale.txt", sep="\n", append=FALSE)
out<-capture.output(Anova(lmf))
cat("Results", out, file="Results/AnovaSummarylmFemale.txt", sep="\n", append=FALSE)


#Linear models for males 
lmm<-lm(diffvar~Course+ResBC+year+MeanLen+logMean+jonvar+drainage, data=FS_trinm) #Linear model for males
vif(lmm)#Checking  male model for inflation factors
Anova(lmm)#Looking at important of each predictor 
summary(lmm)## Summary for linear model for males
out<-capture.output(summary(lmm))
cat("Results", out, file="Results/SummarylmMale.txt", sep="\n", append=FALSE)
out<-capture.output(Anova(lmm))
cat("Results", out, file="Results/AnovaSummarylmMale.txt", sep="\n", append=FALSE)

#Plotting the partial residuals 

#Linear model plot for all individuals ogMeanP
Agg_lmall_logMeanP=visreg(lmallwsex, "logMean", gg=TRUE)+  #Partial residual plot for lm all
  ylab("Difference in Variance")+    #Improved ylab
  xlab("log Mean parasite load")+     #Improved xlab
  theme_classic()+                #Setting theme to classic
  geom_hline(aes(yintercept=0))+  #Setting the hline to 0 to show unity
  theme(text = element_text(size=32))

pdf("Figures/Agg_lmall_logMeanP.pdf")#Saving the figure as a pdf
Agg_lmall_logMeanP
dev.off()

#Linear model plot for all individuals ogMeanP
Agg_lmall_MeanLen=visreg(lmallwsex, "MeanLen", gg=TRUE)+  #Partial residual plot for lm all
  ylab("Difference in Variance")+    #Improved ylab
  xlab("Mean Length (mm)")+     #Improved xlab
  theme_classic()+                #Setting theme to classic
  geom_hline(aes(yintercept=0), size=1)+  #Setting the hline to 0 to show unity
  theme(text = element_text(size=32))

pdf("Figures/Agg_lmall_MeanLen.pdf")#Saving the figure as a pdf
Agg_lmall_MeanLen
dev.off()


#Plot for lmf for logMeanP
Agg_lmf_logMeanParasite=visreg(lmf, "logMean", gg=TRUE, line = list(size=1))+ #Partial residual plot for lm all
  ylab("Difference in Variance")+ #improved ylab
  xlab("Log Mean Parasite Load")+ #Improved xlab
  theme_classic()+#Setting theme to classic
  geom_hline(aes(yintercept=0, size=1.5))+  #Setting the hline to 0 to show unity
  theme(text = element_text(size=32))

pdf("Figures/Agg_lmf_logMeanParasite.pdf")#Saving the figure as a pdf
Agg_lmf_logMeanParasite
dev.off()


#Plot for lmf for logMeanP
Agg_lmf_MeanLen=visreg(lmm, "MeanLen", gg=TRUE)+ #Partial residual plot for lm all
  ylab("Difference in Variance")+ #improved ylab
  xlab("Mean Length (mm)")+ #Improved xlab
  theme_classic()+#Setting theme to classic
  geom_hline(aes(yintercept=0))+  #Setting the hline to 0 to show unity
  theme(text = element_text(size=32))

pdf("Figures/Agg_lmf_MeanLen.pdf")#Saving the figure as a pdf
Agg_lmf_MeanLen
dev.off()


#Plot for lmallwsex for sex
Agg_lmallwsex_sex=visreg(lmallwsex, "sex", gg=TRUE)+ #Partial residual plot for lm all
  ylab("Difference in Variance")+ #improved ylab
  xlab("Sex")+ #Improved xlab
  theme_classic()+#Setting theme to classic
  geom_hline(aes(yintercept=0))+#Setting the hline to 0 to show unity
  theme(text = element_text(size=32))

pdf("Figures/Agg_lmallwsex_sex.pdf")#Saving the figure as a pdf
Agg_lmallwsex_sex
dev.off()

#Plot for lmm for stuff
#Plot for MeanLen for lmm
Agg_lmm_PRjonvar=visreg(lmm, "jonvar", gg=TRUE)+ #Partial residual plot for lm all
  ylab("Difference in Variance")+ #improved ylab
  xlab("Variation in length (Residual Distance)")+ #Improved xlab
  theme_classic()+#Setting theme to classic
  geom_hline(aes(yintercept=0))+#Setting the hline to 0 to show units
  theme(text = element_text(size=24))

pdf("Figures/Agg_lmm_PRjonvar.pdf")#Saving the figure as a pdf
Agg_lmm_PRjonvar
dev.off()
#Plot for MeanLen for lmm
Agg_lmm_PRCourse=visreg(lmm, "Course", gg=TRUE)+ #Partial residual plot for lm all
  ylab("Difference in Variance")+ #improved ylab
  xlab("Predation Regime")+ #Improved xlab
  theme_classic()+#Setting theme to classic
  geom_hline(aes(yintercept=0))+#Setting the hline to 0 to show unity
  theme(text = element_text(size=24))

pdf("Figures/Agg_lmm_PRCourse.pdf")#Saving the figure as a pdf
Agg_lmm_PRCourse
dev.off()

################################ Number of obs vs expected uninfected vs infected analysis
##This part of the analysis got a bit tricky. The models dont exactly love this structure of the data. 
#We tried a regular linear model for just raw differences in zero due to the negatives but also we did that 
#difference plus the absolute minimum to make all values either 0 or positive and the negative binomial model 
#did not like this either. Now we will be moving on to try difference in prevalence given this could potentially work 
#With a beta regression
View(FS_trin)
FS_trin$DiffZero<-as.numeric(FS_trin$DiffZero)
is.na(FS_trin$DiffZero)
FS_trinN<-subset(FS_trin, DiffZero != "NA")


View(FS_trinN)

FS_trinN$a<-min(FS_trinN$DiffZero)*-1
FS_trinN$AdjDZ<-FS_trinN$DiffZero+FS_trinN$a

ZeroModel<-lm(DiffZero~sex+Course+jonvarRS+jonresRS+MeanLenRS+logMeanRS+sex:Course+sex:MeanLenRS, data=FS_trinN)
summary(ZeroModel)

ZeroModel1<-lm(DiffZero~sex+Course+jonvarRS+jonresRS+MeanLenRS+logMeanRS, data=FS_trinN)
summary(ZeroModel1)

ZeroModel2<-glm(AdjDZ~sex+Course+jonvarRS+jonresRS+MeanLenRS+logMeanRS+sex:Course+sex:MeanLenRS, data=FS_trinN)
summary(ZeroModel2)

ZeroModel2<-glmmTMB(AdjDZ~sex+Course+jonvarRS+jonresRS+MeanLenRS+logMeanRS+sex:Course+sex:MeanLenRS+(1|drainage), family=poisson(link="log"),  data=FS_trinN)
summary(ZeroModel2)

ZeroModel3<-glmmTMB(AdjDZ~sex+jonvarRS+jonresRS+MeanLenRS+logMeanRS+sex:Course+sex:MeanLenRS+(1|drainage2), family=nbinom2,  data=FS_trinN)
summary(ZeroModel3)
anova(ZeroModel2,ZeroModel3)


sim_residuals_glmmTMB <-simulateResiduals(ZeroModel1, 1000)  #Quantile residuals
plot(sim_residuals_glmmTMB) 



visreg(ZeroModel)


plot(residuals(ZeroModel, type="pearson")~fitted(ZeroModel)) #Pearson residual plot vs fitted values
hist(residuals(ZeroModel, type="pearson")) #histogram of Pearson residuals
#Variables in the model
plot(residuals(ZeroModel, type="pearson")~FS_trinN$sex) #Residual plot against sex
plot(residuals(ZeroModel, type="pearson")~FS_trinN$Course)#Residual plot against Course
plot(residuals(ZeroModel, type="pearson")~FS_trinN$jonvarRS)#Residual plot against scaled jonvar
plot(residuals(ZeroModel, type="pearson")~FS_trinN$jonresRS)#Residual plot against scaled jonres
plot(residuals(ZeroModel, type="pearson")~FS_trinN$MeanLenRS)#Residual plot against scaled mean length
plot(residuals(ZeroModel, type="pearson")~FS_trinN$logMeanRS)

sim_residuals_glmmTMB <-simulateResiduals(ZeroModel, 1000)  #Quantile residuals
plot(sim_residuals_glmmTMB) 


overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}

overdisp_fun(ZeroModel)


plot(residuals(ZeroModel1, type="pearson")~fitted(ZeroModel1)) #Pearson residual plot vs fitted values
hist(residuals(ZeroModel1, type="pearson")) #histogram of Pearson residuals
#Variables in the model
plot(residuals(ZeroModel1, type="pearson")~FS_trinN$sex) #Residual plot against sex
plot(residuals(ZeroModel1, type="pearson")~FS_trinN$Course)#Residual plot against Course
plot(residuals(ZeroModel1, type="pearson")~FS_trinN$jonvarRS)#Residual plot against scaled jonvar
plot(residuals(ZeroModel1, type="pearson")~FS_trinN$jonresRS)#Residual plot against scaled jonres
plot(residuals(ZeroModel1, type="pearson")~FS_trinN$MeanLenRS)#Residual plot against scaled mean length
plot(residuals(ZeroModel1, type="pearson")~FS_trinN$logMeanRS)

##Analysis for the difference in prevalence. 
#This part of the analysis is using the difference in the observed prevalence of parasites in each population
#versus the expected prevalence of parasites based on the feasible set. We will again use the absolute minimum +0.0001
# plus difference between them to allow all things to be between 0 and 1 to meet the assumptions of the beta distribution

FS_trin$DiffPrev<-FS_trin$ObsPrev - FS_trin$ExpPrev
 summary(FS_trinN$AdjPrev)

hist(FS_trinN$DiffPrev[FS_trinN$sex =="F"], breaks=10)
hist(FS_trinN$DiffPrev[FS_trinN$sex =="M"], breaks = 10)

#PrevM<-glmmTMB(AdjPrev~sex+Course+jonvarRS+jonresRS+logMeanRS+MeanLenRS+sex:MeanLenRS+(1|drainage2),family=beta_family(link="logit") ,data=FS_trinN)
#summary(PrevM)

#sim_residuals_glmmTMB <-simulateResiduals(PrevM, 1000)  #Quantile residuals
#plot(sim_residuals_glmmTMB) 
#testUniformity(sim_residuals_glmmTMB)

#plot(residuals(PrevM, type="pearson")~fitted(PrevM)) #Pearson residual plot vs fitted values
#hist(residuals(PrevM, type="pearson")) #histogram of Pearson residuals
#Variables in the model
#plot(residuals(PrevM, type="pearson")~FS_trinN$sex) #Residual plot against sex
#plot(residuals(PrevM, type="pearson")~FS_trinN$Course)#Residual plot against Course
#plot(residuals(PrevM, type="pearson")~FS_trinN$jonvarRS)#Residual plot against scaled jonvar
#plot(residuals(PrevM, type="pearson")~FS_trinN$jonresRS)#Residual plot against scaled jonres
#plot(residuals(PrevM, type="pearson")~FS_trinN$MeanLenRS)#Residual plot against scaled mean length
#plot(residuals(PrevM, type="pearson")~FS_trinN$logMeanRS)


#LRT 
#PrevM1<-glmmTMB(AdjPrev~sex+jonvarRS+logMeanRS+MeanLenRS+(1|drainage2),family=beta_family(link="logit") ,data=FS_trinN)
#anova(PrevM,PrevM1, test="LRT")

DPM<-glmmTMB(DiffPrev~sex+Course+logMeanRS+MeanLenRS+(1|drainage2),family=gaussian(),FS_trinN,)
summary(DPM)
Anova(DPM)
fixedvar<-varFixed(~logMeanRS)
DPM2<-lme(DiffPrev~sex+Course+logMeanRS+MeanLenRS+sex:MeanLenRS+sex:logMeanRS, random = ~1|drainage2, data=FS_trinN, weight=fixedvar)
summary(DPM2)
Anova(DPM2, type=2)

check_heteroscedasticity(DPM)
AIC(DPM,DPM2)

sim_residuals_DPM <-simulateResiduals(DPM2, 1000)  #Quantile residuals
plot(sim_residuals_DPM) 
testQuantiles(DPM, predictor = FS_trinN$MeanLenRS)
testQuantiles(DPM, predictor = FS_trinN$logMeanRS)
testCategorical(DPM, catPred = FS_trinN$Course)
testCategorical(DPM, catPred = FS_trinN$sex)
testOutliers(sim_residuals_DPM)
testDispersion(sim_residuals_DPM)




t.test(FS_trinN$DiffPrev, mu=0)
t.test(FS_trinN$DiffPrev[FS_trinN$sex=="F"], mu=0)
t.test(FS_trinN$DiffPrev[FS_trinN$sex=="M"], mu=0)
t.test(FS_trinN$DiffPrev[FS_trinN$Course=="Lower"], mu=0)
t.test(FS_trinN$DiffPrev[FS_trinN$Course=="Upper"], mu=0)




plot(residuals(DPM2, type="pearson")~fitted(DPM2)) #Pearson residual plot vs fitted values
hist(residuals(DPM2, type="pearson")) #histogram of Pearson residuals
#Variables in the model
plot(residuals(DPM2, type="pearson")~FS_trinN$sex) #Residual plot against sex
plot(residuals(DPM2, type="pearson")~FS_trinN$Course)#Residual plot against Course
plot(residuals(DPM2, type="pearson")~FS_trinN$jonvarRS)#Residual plot against scaled jonvar
plot(residuals(DPM2, type="pearson")~FS_trinN$jonresRS)#Residual plot against scaled jonres
plot(residuals(DPM2, type="pearson")~FS_trinN$MeanLenRS)#Residual plot against scaled mean length
plot(residuals(DPM2, type="pearson")~FS_trinN$logMeanRS)



visreg(PrevM1)

pairs(~diffvar+sex+MeanW+VarW+Course+
ResBC+jonvar+jonres+MeanLen+logMean+VarLen+
river+year+drainage+AdjPrev, lower.panel=panel.smooth, 
diag.panel=panel.hist, upper.panel=panel.cor, data=FS_trinN)


if (requireNamespace("broom.mixed") && requireNamespace("dotwhisker")) {
  (t1 <- broom.mixed::tidy(DPM2, conf.int = TRUE))
  if (packageVersion("dotwhisker")>"0.4.1") {
    dw <- dwplot(DPM2)
  } else {
    PrevM$coefficients <- TRUE ## hack!
    dw <- dwplot(DPM2,by_2sd=FALSE)
  }
  print(dw+geom_vline(xintercept=0,lty=2))
}



ss <- summary(model2)
## print table; add space,
pxt <- function(x,title) {
  cat(sprintf("{\n\n\\textbf{%s}\n\\ \\\\\\vspace{2pt}\\ \\\\\n",title))
  print(xtable(x), floating=FALSE); cat("\n\n")
  cat("\\ \\\\\\vspace{5pt}\\ \\\\\n")
}

pxt(lme4::formatVC(ss$varcor$cond),"random effects variances")
pxt(coef(ss)$cond,"conditional fixed effects")
pxt(coef(ss)$zi,"conditional zero-inflation effects")

hist(FS_trin$DiffMaxP)




#MaxParasite counts
FS_trinN<-subset(FS_trin, DiffPrev != "NA" & DiffMaxP != "NA")


FS_trinN$Site<-as.factor(FS_trinN$Site)
FS_trinN$drainage<-as.factor(FS_trinN$drainage)
FS_trinN$drainage2<-as.factor(FS_trinN$drainage2)
FS_trinN$river<-as.factor(FS_trinN$river)
FS_trinN$year<-as.factor(FS_trinN$year)
FS_trinN$Course<-as.factor(FS_trinN$Course)
FS_trinN$sex<-as.factor(FS_trinN$sex)




MinmaxP<-min(FS_trinN$DiffMaxP)
hist(FS_trinN$DiffMaxP, breaks=20)
FS_trinN$AdjDiffMaxP<-FS_trinN$DiffMaxP+(MinmaxP*-1)+1
hist(FS_trinN$AdjDiffMaxP, breaks=20)
FS_trinN$AdjDiffMaxP2<-sqrt(FS_trinN$AdjDiffMaxP+1)
FS_trinN1<-subset(FS_trinN, DiffMaxP < 100)

MaxModelM<-lmer(DiffMaxP~sex+Course+jonres+logMean+MeanLen+(1|Site), FS_trinN)
check_model(MaxModelM)
check_distribution(MaxModelM)
check_normality(MaxModelM)

plot(residuals(MaxModelM)~FS_trinN$sex)
summary(MaxModelM)
Anova(MaxModelM)
SexMax<-visreg(MaxModelM, "sex", scale="response", gg=TRUE,  ylab="Difference in largest parasite intensity", xlab="Sex", partial=TRUE) +geom_hline(yintercept=0)+theme_classic()

visreg(MaxModelM, "MeanLenRS", scale="response", partial=TRUE,gg=TRUE, ylab="Difference in largest parasite intensity", xlab="Sex") +geom_hline(yintercept=0)+theme_classic()
t.test(FS_trinN$DiffMaxP[FS_trinN$sex == "M"])

sim_residuals_DPM <-simulateResiduals(MaxModelM, 1000)  #Quantile residuals
plot(sim_residuals_DPM) 

plot(residuals(MaxModelM)~fitted(MaxModelM)) #Pearson residual plot vs fitted values
hist(residuals(MaxModelM)) #histogram of Pearson residuals
#Variables in the model
plot(residuals(MaxModelM)~FS_trinN$sex) #Residual plot against sex
plot(residuals(MaxModelM)~FS_trinN$Course)#Residual plot against Course
plot(residuals(MaxModelM)~FS_trinN$jonvarRS)#Residual plot against scaled jonvar
plot(residuals(MaxModelM)~FS_trinN$jonresRS)#Residual plot against scaled jonres
plot(residuals(MaxModelM)~FS_trinN$MeanLenRS)#Residual plot against scaled mean length
plot(residuals(MaxModelM)~FS_trinN$logMeanRS)




#Some stats to clean up our data

summary(aov(lm(diffvar~year, FS_trin))) #diffvar doesnt vary temporally
summary(aov(lm(DiffPrev~year, FS_trinN))) #diffvar doesnt vary temporally
summary(aov(lm(DiffMaxP~year, FS_trinN)))


corrcsv<-read.csv("FeasCorr")

ggplot(FeasCorr, aes(x=FemaleDiff, y=MaleDiff))+geom_point()+geom_abline( slope=1)
library(ggplot2)

cor(FeasCorr$FemaleDiff,FeasCorr$MaleDiff)
