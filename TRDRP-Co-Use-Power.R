#TRDRP Co-Use Grant Application
#PI: Spencer Bujarski
#Power Calculation via simulation

library(SpPack)
library(lme4)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gmodels) #CrossTable Function
library(Hmisc) #rcorr

#Want to see if Alcohol/Marijuana use increases chances of cigarette smoking
#all variables are binary
#plan is 2-weeks of daily data collection per subject
#EMA study with lots of observations, but for simulation just treat each day as a single potential drug use event

#SubjectSim function
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

#run full simulation to test SubjectSim
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


#Simulate multiple subjects----
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



#Simulate with heterogeneous parameters at subject level----

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


#Run full power calculation via simulation with nested for loops----
#number of simulated datasets set to 1000 per parameter combination
Nsims <- 1000

#Start with power simulation of Odds Ratio = 2
#Multiple variable parameters for power simulation
#NSubs - number of subjects between 50 and 300 by 20
#baseCig - baseline cigarette smoking rate on non-drinking days between .2 and .8 by .2
#baseAlc - baseline drinking rate on non-drinking days between .2 and .8 by .2
#AlcOR - Odds ratio of Alcohol effect set to 2
#Days - days of observation set to 14
#Returns PowerSim.OR2 dataset with outcome variables of power.05 and power.01 based on simulated power with alpha = 0.05 and 0.01

PowerSim.OR2 <- data.frame(expand.grid(NSubs = seq(50, 300, 20),
                                       baseCig = seq(.20, .80, .20),
                                       AlcRate = seq(.20, .80, .20),
                                       AlcOR = 2,
                                       Days = 14,
                                       Power.05 = NA,
                                       Power.01 = NA))

#Outer for loop to walk through power calculation parameters in PowerSim.OR2
for(p in 1:dim(PowerSim.OR2)[1]) { 
  
  #vector of pvalues to store
  pvalues <- rep(NA, Nsims)
  
  #Inner for loop to generate NSims simulated datasets, perform analysis and store pvalues
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
    
    #This older approach might repeat p-values if the new model doesn't converge for some reason. 
    # logistic.MLM <- glmer(Cig ~ (1 | Subject) + Alc,
    #                       data=FullData,
    #                       family=binomial,
    #                       control = glmerControl(optimizer = "bobyqa"))
    
    #extract Alc p-value for power calculation later
    pvalues[n] <- summary(glmer(Cig ~ (1 | Subject) + Alc, data=FullData,  
                                family=binomial, control = glmerControl(optimizer = "bobyqa")))$coefficients["Alc","Pr(>|z|)"]
  }
  
  #calculate power based on simulated datasets and analysis
  PowerSim.OR2$Power.05[p] <- sum(pvalues < 0.05)/Nsims #sum(!is.na(pvalues))
  PowerSim.OR2$Power.01[p] <- sum(pvalues < 0.01)/Nsims
}
  
#plotting OR2 power analysis
PowerSim.OR2$AlcRate.Str <- paste("Drinking Rate: ", PowerSim.OR2$AlcRate)
  
colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=4))
PowerSim.OR2.plot.05 <- ggplot(PowerSim.OR2, aes(x = NSubs, y=Power.05, colour = as.factor(baseCig))) +
  geom_line(size = 2) +
  facet_wrap(~ AlcRate.Str) +
  scale_colour_manual("Baseline Smoking Rate", values=colorscale) + 
  scale_x_continuous("Sample Size") + 
  scale_y_continuous("Power (at alpha = 0.05)") +
  ggtitle("Power from simulation with Odds Ratio = 2, Alpha = 0.05\nLevel 2 heterogeneous") +
  DotRTheme(legend.position = "right", title.size = 16)
PowerSim.OR2.plot.05
ggsave(PowerSim.OR2.plot.05, filename="PowerSim.OR2.plot.05.png", width = 8, height=6, dpi=250)

colorscale <- scales::seq_gradient_pal("lightblue", "navyblue", "Lab")(seq(0,1,length.out=4))
PowerSim.OR2.plot.01 <- ggplot(PowerSim.OR2, aes(x = NSubs, y=Power.01, colour = as.factor(baseCig))) +
  geom_line(size = 2) +
  facet_wrap(~ AlcRate.Str) +
  scale_colour_manual("Baseline Smoking Rate", values=colorscale) + 
  scale_x_continuous("Sample Size") + 
  scale_y_continuous("Power (at alpha = 0.01)") +
  ggtitle("Power from simulation with Odds Ratio = 2, Alpha = 0.01\nLevel 2 heterogeneous") +
  DotRTheme(legend.position = "right", title.size = 16)
PowerSim.OR2.plot.01
ggsave(PowerSim.OR2.plot.01, filename="PowerSim.OR2.plot.01.png", width = 8, height=6, dpi=250)



#Make a function to compile hopefully to speed up
TRDRP.PowerSim <- function(Nsims, NSubs, baseCig, AlcRate, AlcOR, Days){
  #create dataframe of parameters to test
  PowerSim <- expand.grid(NSubs = NSubs,
                          baseCig = baseCig,
                          AlcRate = AlcRate,
                          AlcOR = AlcOR,
                          Days = Days,
                          Power.05 = NA,
                          Power.01 = NA)
  
  #Outer for loop to walk through power calculation parameters in PowerSim
  for(p in 1:dim(PowerSim)[1]) {

    #vector of pvalues to store
    pvalues <- rep(NA, Nsims)

    #Inner for loop to generate NSims simulated datasets, perform analysis and store pvalues
    for (n in 1:Nsims){ #run Nsims data simulations for power calculations
      #print where sim is
      print(noquote(paste("PowerSim", p, " of ", dim(PowerSim)[1], " Sim number", n, " of ", Nsims)))

      #fixed parameters
      NSubs.t <- PowerSim$NSubs[p] #pull NSubs.t from PowerSim dataframe
      Days.t <- PowerSim$Days[p] # same for Days.t

      #Simulate randomly variable subject parameters -- baseCig.t, AlcRate.t, AlcOR.t
      parameters <- data.frame(baseCig.t = rnorm(NSubs.t, PowerSim$baseCig[p], .2), #simulate base Cig smoking rate based on PowerSim. SD fixed at .2
                               AlcRate.t = rnorm(NSubs.t, PowerSim$AlcRate[p], .2), #same for AlcRate.t
                               AlcOR.t = rlnorm(n=NSubs.t, meanlog=log(PowerSim$AlcOR[p]), sdlog=.4)) #log normal distribution of subject Odds Ratio (mean OR from PowerSim, sd=.4)

      #Ensure no negative rates or >=100% rates
      #will slightly bias results to the mean, but whatever
      while(min(parameters$baseCig.t) <= 0 || max(parameters$baseCig.t) >= 1){
        for(i in 1:NSubs.t){
          if(parameters$baseCig.t[i] <= 0 || parameters$baseCig.t[i] >= 1){parameters$baseCig.t[i] <- rnorm(1, PowerSim$baseCig[p], .2)}
        }
      }
      while(min(parameters$AlcRate.t) <= 0 || max(parameters$AlcRate.t) >= 1){
        for(i in 1:NSubs.t){
          if(parameters$AlcRate.t[i] <= 0 || parameters$AlcRate.t[i] >= 1){parameters$AlcRate.t[i] <- rnorm(1, PowerSim$AlcRate[p], .2)}
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

      #extract Alc p-value for power calculation later
      pvalues[n] <- summary(glmer(Cig ~ (1 | Subject) + Alc, data=FullData,
                                  family=binomial, control = glmerControl(optimizer = "bobyqa")))$coefficients["Alc","Pr(>|z|)"]
    }

    #calculate power based on simulated datasets and analysis
    PowerSim$Power.05[p] <- sum(pvalues < 0.05)/sum(!is.na(pvalues))
    PowerSim$Power.01[p] <- sum(pvalues < 0.01)/sum(!is.na(pvalues))
  }
  
  return(PowerSim)
}


TRDRP.PowerSim(Nsims=10, NSubs=seq(50, 100, 25), baseCig=c(.2, .5), AlcRate=c(.2, .5), AlcOR=2, Days=14)

#test speed
system.time(TRDRP.PowerSim(Nsims=10, NSubs=seq(50, 100, 25), baseCig=c(.2, .5), AlcRate=c(.2, .5), AlcOR=2, Days=14))
#user  system elapsed 
#30.95    0.11   33.12 

require(compiler)
cmpfun(TRDRP.PowerSim)
system.time(TRDRP.PowerSim(Nsims=10, NSubs=seq(50, 100, 25), baseCig=c(.2, .5), AlcRate=c(.2, .5), AlcOR=2, Days=14))
#   user  system elapsed 
#28.48    0.06   29.61 
#Not that much speedup

cmpfun(SubjectSim)
system.time(TRDRP.PowerSim(Nsims=10, NSubs=seq(50, 100, 25), baseCig=c(.2, .5), AlcRate=c(.2, .5), AlcOR=2, Days=14))
#   user  system elapsed 
#31.50    0.06   31.59
#No benefit at all of compiler

PowerSim.OR15 <- TRDRP.PowerSim(Nsims=1000, NSubs=seq(50, 300, 25), baseCig=c(.2, .5), AlcRate=c(.2, .5), AlcOR=1.5, Days=14)

