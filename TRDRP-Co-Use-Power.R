#TRDRP Co-Use Grant Application
#PI: Spencer Bujarski
#Power Calculation via simulation

library(SpPack)
library(lme4)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gmodels) #CrossTable Function
library(Hmisc) #rcorr

#Want to see if Alcohol/Marijuana use increases chances of cigarette smoking
#all variables are binary
#plan is 2-weeks of daily data collection per subject
#EMA study with lots of observations, but for simulation just treat each day as a single potential drug use event

#Function to generate single subject data
#baseCig = baseline smoking probability within subject
#AlcRate = baseline drinking probability within subject
#AlcOR = Odds ratio for alcohol effect within this subject
#Days = Days of observation (e.g. n potential drug using events)
SubjectSim <- function(baseCig, AlcRate, AlcOR, Days){
  
  Alc <- rbinom(Days, 1, AlcRate) #simulate 14 days of possible drinking based on AlcRate
  
  ObaseCig <- baseCig / (1 - baseCig) #baseline odds of cig
  AlcCig <- (AlcOR*ObaseCig)/(1+AlcOR*ObaseCig)
  
  #printing for testing
  #print(paste("ObaseCig",ObaseCig))
  print(paste("AlcCig",AlcCig))
  
  Cig <- ifelse(Alc==0, rbinom(sum(Alc==0), 1, baseCig), #simulate cigarette smoking with no alcohol
                rbinom(sum(Alc==1), 1, AlcCig)) #simulate cigarette smoking after drinking
    
  DaysData <- data.frame(Alc=Alc, Cig=Cig)
  return(DaysData)
}

#SubjectSim Testing
SubjectSim(baseCig=.1, AlcRate = .3, AlcOR = 9, Days = 14)

#run full simulation
sim=1000
NonDrinkCigMeans <- rep(NA, sim)
DrinkCigMeans <- rep(NA, sim)
for(i in 1:sim){
  SubData <- SubjectSim(baseCig=.1, AlcRate = .4, AlcOR = 5, Days = 14)
  
  CigMeans <- SubData %>% group_by(Alc) %>% summarise(CigMean=mean(Cig))
  
  NonDrinkCigMeans[i] <- CigMeans$CigMean[1]
  DrinkCigMeans[i] <- CigMeans$CigMean[2]
}

#seems to work great!
SpHist(NonDrinkCigMeans)
SpHist(DrinkCigMeans)


#simulate multiple subjects
NSubs.t = 100 #test parameters
baseCig.t = .50
AlcRate.t = .50
AlcOR.t = 2
Days.t = 14

#first try with homogenous population odds ratio
FullData <- cbind(Subject=1, SubjectSim(baseCig=baseCig.t, AlcRate = AlcRate.t, AlcOR = AlcOR.t, Days = Days.t))
for (i in 2:NSubs.t){
  SubjectData <- cbind(Subject=i, SubjectSim(baseCig=baseCig.t, AlcRate = AlcRate.t, AlcOR = AlcOR.t, Days = Days.t))
  FullData <- rbind(FullData, SubjectData)
}

FullData %>% group_by(Subject,Alc) %>% summarise(CigMean=mean(Cig))

SpHist(subset(FullData %>% group_by(Subject,Alc) %>% summarise(CigMean=mean(Cig)), Alc == 0)$CigMean)
SpHist(subset(FullData %>% group_by(Subject,Alc) %>% summarise(CigMean=mean(Cig)), Alc == 1)$CigMean)

#run a sample model with simulated subjects
#Seems to work great using lme4:glmer
#https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
logistic.MLM <- glmer(Cig ~ (1 | Subject) + Alc,
                    data=FullData,
                    family=binomial,
                    control = glmerControl(optimizer = "bobyqa"))
summary(logistic.MLM)



#Simulate with subject specific values

#universal parameters
NSubs.t = 100 #test parameters
Days.t = 14

#subject specific parameters
parameters <- data.frame(baseCig.t = rnorm(NSubs.t, .5, .2), #simulate base Cig smoking rate from N(.5, .2)
                         AlcRate.t = rnorm(NSubs.t, .5, .2), #simulate base Cig smoking rate from N(.5, .2))
                         AlcOR.t = rlnorm(n=NSubs.t, meanlog=log(2), sdlog=.4)) #log normal distribution of subject Odds Ratio (mean OR = 2, sd=.4)
SpDesc(parameters)
#check for no negative rates
while(min(parameters$baseCig.t) <=0){
  parameters$baseCig.t <- ifelse(parameters$baseCig.t<=0, rnorm(NSubs.t, .5, .2), parameters$baseCig.t)
}

while(min(parameters$AlcRate.t) <=0){
  parameters$AlcRate.t <- ifelse(parameters$AlcRate.t<=0, rnorm(NSubs.t, .5, .2), parameters$AlcRate.t)
}
SpDesc(parameters)
SpHist(parameters, variable = "baseCig.t")
SpHist(parameters, variable = "AlcRate.t")
SpHist(parameters, variable = "AlcOR.t")

#simulate Nsubs.t subject data
FullData <- cbind(Subject=1, SubjectSim(baseCig=parameters$baseCig.t[1], AlcRate = parameters$AlcRate.t[1], 
                                        AlcOR = parameters$AlcOR.t[1], Days = Days.t))
for (i in 2:NSubs.t){
  SubjectData <- cbind(Subject=i, SubjectSim(baseCig=parameters$baseCig.t[i], AlcRate = parameters$AlcRate.t[i], 
                                             AlcOR = parameters$AlcOR.t[i], Days = Days.t))
  FullData <- rbind(FullData, SubjectData)
}

FullData %>% group_by(Subject,Alc) %>% summarise(CigMean=mean(Cig))

logistic.MLM <- glmer(Cig ~ (1 | Subject) + Alc,
                      data=FullData,
                      family=binomial,
                      control = glmerControl(optimizer = "bobyqa"))
summary(logistic.MLM)

#extract Alc p-value for power calculation later
summary(logistic.MLM)$coefficients["Alc","Pr(>|z|)"]


#Now that simulation works to make a MLM dataset, need to write some for loops to run simulations for power analysis





