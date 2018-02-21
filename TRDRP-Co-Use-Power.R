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
  #print(paste("AlcCig",AlcCig))
  
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

#Simulate randomly variable subject parameters -- baseCig.t, AlcRate.t, AlcOR.t
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


#Now that simulation works to make a nested dataset, need to write some for loops to run simulations for power analysis

PowerSim.OR2 <- data.frame(expand.grid(NSubs = seq(50, 300, 20),
                                       baseCig = seq(.20, .80, .20),
                                       AlcRate = seq(.20, .80, .20),
                                       AlcOR = 2,
                                       Days = 14,
                                       Power.05 = NA,
                                       Power.01 = NA))
dim(PowerSim.OR2)

Nsims <- 1000
for(p in 1:dim(PowerSim.OR2)[1]) { #walk through power calculations varying NSubs, baseCig, and AlcRate
  
  #vector of pvalues to store
  pvalues <- rep(NA, Nsims)
  for (n in 1:Nsims){ #run Nsims data simulations for power calculations
    #print where sim is
    print(noquote(paste("PowerSim", p, " of ", dim(PowerSim.OR2)[1], " Sim number", n, " of ", Nsims)))
   
  
    #fixed parameters
    NSubs.t <- PowerSim.OR2$NSubs[p] #pull NSubs.t from PowerSim.OR2 dataframe
    Days.t <- PowerSim.OR2$Days[p] # same for Days.t
    
    #Simulate randomly variable subject parameters -- baseCig.t, AlcRate.t, AlcOR.t
    parameters <- data.frame(baseCig.t = rnorm(NSubs.t, PowerSim.OR2$baseCig[p], .2), #simulate base Cig smoking rate based on PowerSim.OR2. SD fixed at .2
                             AlcRate.t = rnorm(NSubs.t, PowerSim.OR2$AlcRate[p], .2), #same for AlcRate.t
                             AlcOR.t = rlnorm(n=NSubs.t, meanlog=log(PowerSim.OR2$AlcOR[p]), sdlog=.4)) #log normal distribution of subject Odds Ratio (mean OR from PowerSim.OR2, sd=.4)
    
    #Ensure no negative rates or >=100% rates
    #will slightly bias results to the mean, but whatever
    while(min(parameters$baseCig.t) <= 0 || max(parameters$baseCig.t) >= 1){
      for(i in 1:NSubs.t){
        if(parameters$baseCig.t[i] <= 0 || parameters$baseCig.t[i] >= 1){parameters$baseCig.t[i] <- rnorm(1, PowerSim.OR2$baseCig[p], .2)}
      }
    }
    while(min(parameters$AlcRate.t) <= 0 || max(parameters$AlcRate.t) >= 1){
      for(i in 1:NSubs.t){
        if(parameters$AlcRate.t[i] <= 0 || parameters$AlcRate.t[i] >= 1){parameters$AlcRate.t[i] <- rnorm(1, PowerSim.OR2$AlcRate[p], .2)}
      }
    }
    
    #simulate full nested dataset
    FullData <- cbind(Subject=1, SubjectSim(baseCig=parameters$baseCig.t[1], AlcRate = parameters$AlcRate.t[1], 
                                            AlcOR = parameters$AlcOR.t[1], Days = Days.t))
    for (i in 2:NSubs.t){
      SubjectData <- cbind(Subject=i, SubjectSim(baseCig=parameters$baseCig.t[i], AlcRate = parameters$AlcRate.t[i], 
                                                 AlcOR = parameters$AlcOR.t[i], Days = Days.t))
      FullData <- rbind(FullData, SubjectData)
    }
    
    logistic.MLM <- glmer(Cig ~ (1 | Subject) + Alc,
                          data=FullData,
                          family=binomial,
                          control = glmerControl(optimizer = "bobyqa"))
    
    #extract Alc p-value for power calculation later
    pvalues[n] <- summary(logistic.MLM)$coefficients["Alc","Pr(>|z|)"]
  }
  
  PowerSim.OR2$Power.05[p] <- sum(pvalues < 0.05)/Nsims
  PowerSim.OR2$Power.01[p] <- sum(pvalues < 0.01)/Nsims
}
  

