###############################
##Reads in Required Libraries##
################################
library(multcomp)
library(ggplot2)
library(lmtest)
library(gstat)
library(car)
library(MASS)
library(nlme)
source("C:/Users/Rental/Documents/Stat469/sourceCode/stdres.gls.R")
source("C:/Users/Rental/Documents/Stat469/sourceCode/predictgls.R")
##################
#Reads in Data####
##################
housingData = read.csv(file = "C:/Users/Rental/Documents/Stat469/HousingPrices.csv", header = TRUE, 
              dec = ".")
#Turns Categorical Variables into Factors
housingData$House.Style = as.factor(housingData$House.Style)
housingData$Central.Air = as.factor(housingData$Central.Air)
#Separated data into one containing Na and one with no NA
housingNaData =  housingData[rowSums(is.na(housingData)) > 0,]
noNaData =  na.omit(housingData)
#######
#EDA###
#######
#Central Air Boxplot
unique(noNaData$Central.Air)
YCentral = noNaData[noNaData$Central.Air== "Y",]
NCentral = noNaData[noNaData$Central.Air == "N",]
airNames = c("Yes","No")
boxplot(YCentral$Price,NCentral$Price,ylab = "Appraisal Price",
        main = "Central Air and Price box plot", names = airNames)
#Scatterplot of Above ground living area and Price
ggplot(data= noNaData,mapping=aes(x= Gr.Liv.Area, y= Price)) + geom_point() 
summary(noNaData$Price)
##Garage Car Capacity and Price
#Note that homes with garages that can fit 4 cars were excluded since there was only one home that
#fit this category.
unique(noNaData$Garage.Cars)
garage0  = noNaData[noNaData$Garage.Cars== 0,]
garage1  = noNaData[noNaData$Garage.Cars== 1,]
garage2  = noNaData[noNaData$Garage.Cars== 2,]
garage3  = noNaData[noNaData$Garage.Cars== 3,]
garage4  = noNaData[noNaData$Garage.Cars== 4,]
carNames = c("0","1","2","3")
boxplot(garage0$Price,garage1$Price,garage2$Price,garage3$Price,ylab = "Appraisal Price",
        main = "Garage Car Capacity and Price box plot", names = carNames)
#########################
#Fits appropriate model#
#########################
#Fits MLR model
housinglm =  lm(formula = Price ~ ., data = noNaData)
summary(housinglm)
###############################
##Checking For equal Variance##
###############################
#Using the BP test, I obtained a p-value of less than 2.2e-16. This means that we reject the Null hypothesis
#which means that the assumption of equal variance is not met.
yPred  = fitted(housinglm)
plot(yPred,stdres(housinglm), main= "Fitted values vs Standardized residuals plot",
     xlab = "Fitted values", ylab = "Standardized Residuals", pch = 19)
abline(h = 0, col = 'red', lwd = 2)
bptest(housinglm)
#############################
#Checking for correlation####
#############################
#Variogram:  Variogram is non-linear. This suggests that there is spatial Correlation
myVariogram = variogram(object = Price ~ Gr.Liv.Area + House.Style + Year.Remod.Add + 
               Central.Air + Full.Bath + Half.Bath + Bedroom.AbvGr + Garage.Cars, 
               locations = ~Lon+Lat, data = noNaData)
plot(myVariogram,main = "Variogram of the residuals")
##############################
##Model Fitting###############
##############################
housingExp = gls(model = Price ~ . -Lon - Lat, 
                 data = noNaData,weights = varExp(form = ~ Gr.Liv.Area + Year.Remod.Add + Garage.Cars + Bedroom.AbvGr)
                 ,correlation = corExp(form = ~ Lon + Lat, nugget = TRUE) ,method = "ML")
housingSph = gls(model = Price ~ . -Lon -Lat , 
                  data = noNaData,weights = varExp(form = ~ Gr.Liv.Area + Year.Remod.Add + Garage.Cars + Bedroom.AbvGr),
                  correlation = corSpher(form = ~Lon+Lat, nugget=TRUE) ,  method = "ML")
housingGauss = gls(model = Price ~ -Lon -Lat, 
                  data = noNaData,weights = varExp(form = ~ Gr.Liv.Area + Year.Remod.Add + Garage.Cars + Bedroom.AbvGr),
                  correlation = corGaus(form = ~Lon+Lat, nugget=TRUE) ,  method = "ML")
###########################
##Picking the best Model###
###########################
# Exponential is the best model
AIC(housingExp)
AIC(housingSph)
AIC(housingGauss)
############################
#Model Validation###########
############################
#Linearity Assumption
avPlots(housinglm)
#Independece Assumption
decResids = stdres.gls(housingExp)
residDF = data.frame(Lon = noNaData$Lon,Lat = noNaData$Lat, decorrResid = decResids)
residVariogram  =  variogram(object = decorrResid ~1,locations = ~Lon+Lat, data = residDF)
plot(residVariogram)
#Normality Assumption
hist(decResids)
ks.test(decResids,"pnorm")
#Equal Variance
decResids = stdres.gls(housingExp)
fittedVals = fitted(housingSph)
plot(fittedVals,decResids, main= "Fitted values vs Decorrelated residuals plot",
     xlab = "Fitted values", ylab = "Standardized Residuals", pch=19)
abline(h = 0, col = 'red')
###########################
##Cross-validation########
##########################
#Cross validation
nCV = 1
nTest = 45  #Number of observations in a test set
rpmse = rep(x = NA, times = nCV)
bias = rep(x = NA, times = nCV)
wid = rep(x = NA, times = nCV)
cvg = rep(x = NA, times = nCV)
pb = txtProgressBar(min = 0, max = nCV, style = 3)
for(i in 1:nCV) {
  #Selects tests observations
  testObs = sample(x = 1:nrow(noNaData), size = nTest)
  #Split into training and test data sets
  testData = noNaData[testObs,]
  trainData = noNaData[-testObs,]
  #Fits a linear model using the training data
  glsTrain = gls(model = Price ~ . -Lon - Lat, 
                 data = trainData,weights = varExp(form = ~ Gr.Liv.Area),
                 correlation = corExp(form = ~ Lon + Lat, nugget = TRUE) ,method = "ML")
  #Generates prediction for the test set
  yPredgls = predictgls(glsobj = glsTrain, newdframe = testData, level = 0.95)
  #calculates bias
  bias[i] = mean(yPredgls[,'Prediction'] - testData[['Price']])
  #Calculates RPMSE
  rpmse[i] = (testData[['Price']] - yPredgls[,'Prediction'])^2 %>% mean() %>% sqrt()
  #Calculates coverage
  cvg[i] = ((testData[['Price']] > yPredgls[,'lwr']) & (testData[['Price']] < yPredgls[,'upr'])) %>% mean()
  #Calculates Width
  wid[i] = (yPredgls[,'upr'] - yPredgls[,'lwr']) %>% mean()
  ## Update the progress bar
  setTxtProgressBar(pb, i)
}
close(pb)
#Diagnostic results
#Bias of prediction intervals: 
mean(bias)
#RPMSE of prediction intervals: 
mean(rpmse)
#Coverage of prediction intervals: 
hist(cvg)
mean(cvg)
#Width of prediction intervals: gls mdodel had a narrower confidence interval at 1.92 as opposed to
#the 5.39 for the lm model
mean(wid)
#########################
#Research Questions 1#####
#########################
#Rsquared
SSE = (fitted(housingExp)-noNaData$Price)^2 %>% sum()
SST = (noNaData$Price - mean(noNaData$Price))^2 %>% sum()
1 - (SSE/SST)
#Pseudo-Rsquared
cor(fitted(housingExp),noNaData$Price)^2
#########################
#Research Questions 2#####
#########################
summary(housingExp)
#########################
#Research Questions 3#####
#########################
intervals(housingExp, level = 0.95)[3]
coef(housingExp$modelStruct, unconstrained = FALSE)
#########################
#Research Questions 4#####
#########################
housingNaData =  housingData[rowSums(is.na(housingData)) > 0,]
yPredicted = predictgls(glsobj = housingExp, newdframe = housingNaData,level = 0.95)
housingNaData$Price = yPredicted$Prediction
head(housingNaData)



