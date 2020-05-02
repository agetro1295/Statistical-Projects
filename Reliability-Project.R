#####################
##Required Library##
#####################
library(R2jags)
library(ggplot2)
#####################
##Required Functions#
#####################
GoF_Test = function(fitted_quantiles) {
  n  = length(fitted_quantiles)
  K  = round((n)^(0.4))
  mK = table(cut(fitted_quantiles,(0:K)/K))
  np = n/K
  RB = sum(((mK-np)^2)/np)
  return(1-pchisq(RB,K-1))
}
######################
#Reads in Data#######
####################
cupsData = read.table("C:/Users/rpo22/Documents/466Files/data/cupsData.txt",
                      header = TRUE,sep = ",")
####################
#Data preperation###
####################
cupsData$Sealant_type = as.character(cupsData$Sealant_type)
cupsData$Seal[cupsData$Sealant_type == "F"] = "Flex-seal"
cupsData$Seal[cupsData$Sealant_type == "R"] = "Ruste-oleum"
cupsData$Sealant_type[cupsData$Sealant_type == "F"] = "0"
cupsData$Sealant_type[cupsData$Sealant_type == "R"] = "1"
cupsData$Sealant_type = as.integer(cupsData$Sealant_type)
###########
###EDA#####
###########
failures = cupsData[cupsData$Result == 0 ,]
successes = cupsData[cupsData$Result == 1 ,]
#Histogram of failure by Hole size
hist(failures$Hole_Size, main = "Number of Failures by Hole Size", xlab = "Hole Size")
#Histogram of failures
ggplot(failures, aes(x = Sealant_type,color = factor(Seal))) +
  geom_histogram(binwidth = 0.3) + 
  scale_color_discrete(name = "Sealant Type") + 
  ggtitle("Number of Failures by Sealant Type")
#Setting up the variables we will use in JAGS
n =  length(cupsData[,1])
yFailure = cupsData[,3]
sealantType = cupsData[,4]
holeSize = cupsData[,2]
# Defines the logit and inverse logit functions
logit = function(x) {log(x/(1-x))}
ilogit = function(x) {1/(1+exp(-x))}
#Setting up the Model
cupModel = "model {
  for(i in 1:n){
    yFailure[i] ~ dbin(pi[i],1)
    logit(pi[i]) = beta[1] + beta[2]*holeSize[i] + beta[3]*sealantType[i] 
  }
  beta[1] ~ dnorm(0,1/1000)
  beta[2] ~ dnorm(0,1/1000)
  beta[3] ~ dnorm(0,1/1000)
}
"
#Performs MCMC Simulation
cupSim = jags(
  data=c('yFailure','sealantType','holeSize','n'),
  parameters.to.save=c('beta'),
  model.file=textConnection(cupModel),
  n.iter=21000,
  n.burnin=1000,
  n.chains=5,
  n.thin=1
) 
############################
##Checking For Convergence#
###########################
head(cupSim$BUGSoutput$sims.matrix)
intercept  = cupSim$BUGSoutput$sims.matrix[,1]
holeSizeEffect  = cupSim$BUGSoutput$sims.matrix[,2]
sealantTypeEffect = cupSim$BUGSoutput$sims.matrix[,3]
#Checks mixing: Mixing looks good
par(mfrow=c(3,1))
plot(intercept,type="l")
plot(holeSizeEffect ,type="l")
plot(sealantTypeEffect,type="l")
par(mfrow=c(1,1))
#Effective sample size: Effective sample size is the same as obtained sample
effectiveSize(cupSim$BUGSoutput$sims.matrix[,1])
effectiveSize(cupSim$BUGSoutput$sims.matrix[,2])
effectiveSize(cupSim$BUGSoutput$sims.matrix[,3])
#Gelman Diagnostic Test:
gelman.diag(cupSim$BUGSoutput)
########################
##Goodness of fit test##
#######################
GoF = matrix(NA,ncol=length(yFailure),nrow=length(intercept))
for (i in 1:length(intercept)) {
  GoF[i,] = runif(length(yFailure),
                   pbinom(yFailure-1,1,ilogit(intercept[i] + holeSizeEffect[i]*holeSize 
                    + sealantTypeEffect[i]*sealantType)),
                  pbinom(yFailure,1,ilogit(intercept[i] + holeSizeEffect[i]*holeSize 
                    + sealantTypeEffect[i]*sealantType))
  )
}
# Calculating the p-values for each posterior model
GoF_Summary = apply(GoF,1,GoF_Test)
# Histogram of posterior model p-values
hist(GoF_Summary,xlim=c(0,1))
# Percent of posterior models with p-value less than 0.05
mean(GoF_Summary < 0.05)
####
#Research Question 1
quantile(sealantTypeEffect ,c(0.025,0.975))
hist(sealantTypeEffect, main = "Beta parameter Distribution for Sealant")
########
#PPD###
#######
piFlexSeal = ilogit(intercept + holeSizeEffect*1 + sealantTypeEffect*0)
PPDFlex  =  rbinom(length(piFlexSeal),1,piFlexSeal)
piRust = ilogit(intercept + holeSizeEffect*1 + sealantTypeEffect*1)
PPDRust = rbinom(length(piRust),1,piRust)
sum(PPDRust)/length(PPDRust)
sum(PPDFlex)/length(PPDFlex)
#Research Question 2
quantile(holeSizeEffect,c(0.025,0.975))
hist(holeSizeEffect , main = "Beta parameter Distribution for Hole Size")
#Interpretation of B1
mean(ilogit(intercept + holeSizeEffect*0.5)) -
mean(ilogit(intercept + holeSizeEffect*1.0))
