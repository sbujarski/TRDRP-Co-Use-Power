#TRDRP Co-Use Grant Application
#PI: Spencer Bujarski
#Power Calculation via simulation

library(SpPack)
library(multilevel)
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
  
  print(paste("ObaseCig",ObaseCig))
  print(paste("AlcCig",AlcCig))
  
  Cig <- ifelse(Alc==0, rbinom(1, 1, baseCig), #simulate cigarette smoking with no alcohol
                rbinom(1, 1, AlcCig)) #simulate cigarette smoking after drinking
    
  DaysData <- data.frame(Alc=Alc, Cig=Cig)
  return(DaysData)
}

SubjectSim(baseCig=.1, AlcRate = .4, AlcOR = 9, Days = 14)

for(i in 1:100){
  SubData <- SubjectSim(baseCig=.1, AlcRate = .4, AlcOR = 9, Days = 14)
  
  SubData %>% group_by(as.factor(Alc)) %>% summarise(CigMean=mean(Cig))
  
  CigMeans <-
}



