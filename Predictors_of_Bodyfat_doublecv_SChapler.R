#######STEP 1 DATA REVIEW AND RELATIONSHIP ANALYSIS###################################
######################################################################################

###1.a  Read in file and review data structure.  Look for missing values#############
#read in file and summarize
ais=read.csv("ais.csv")
summary(ais)

#establish dimensions of dataset - 202 rows, 13 variables
dim(ais)
#List names of variables
names(ais)
summary(ais[,1])

######1.b Review variable distribution, normality, equal variance################### 
#Review for normality and outliers

#plot all variables review for relationship with Bfat
par(mfrow=c(1,1))
for (i in (1:13)) {
  plot(ais[,i], ais$Bfat, xlab=paste(names(ais)[i]), ylab="Bodyfat", las=2 )
}


#plot all variables review for influential outliers.  variance across variables appears similar
for (i in 4:13) {
  boxplot(ais[,i], ylab=paste(names(ais)[i]), las=2 )
}

#review for gaussian distribution - see impact of outliers.  SSF right skewed.
for (i in 4:13) {
  hist(ais[,i], main = paste(names(ais)[i]))
}

#no distinct groups - uneven distribution across factors.  unequal variance
plot(ais$Sport)

corrgram(ais, main="Pearson Correlation Between AIS Variables")
       

####1.c Check for signficant relationships, impact of outliers, and collinearity amongst variables##################################
#Full Data fit  R2=.9864, AdjR2=.9857 (non logSSF)
#Full Data fit R2=.9887 AdjR2=0.988
lmfit=lm(Bfat~., data=ais) 
lmfit
summary(lmfit) #clearly SSF, LBM, Wt,Sex have significant relationships with body fat

####variable importance review########  LBM, Wt, SSF, Sex.  Sport levels not strong with over 50% week
#install.packages('caret')
library(caret) 
varImp(lmfit)

 

###see high leverage points in data 187, 195, 200
par(mfrow=c(1,2))
plot(lmfit)
par(mfrow=c(1,1))
###check for collinearity#####
library(car)
vif(lmfit) #high correlation Ht, wt, LBM, BMI, SSF and Hc/Hg

#SSF & Bfat - strong linear relationship
#BMI & Wt, BMI&LBM, LBM&wt, ht& wt, ht & LBM - all appear to have linear relationships with each other and potentially BFat
#HG & HC, HG&RCC, HC & RCC - linear relationships with each other, but not Bfat
par(mfrow=c(1,1))
pairs(ais[,-1])
#better visual of relationships with corrgram with log transformation
#install.packages('corrgram')

library(corrgram)
corrgram(ais, main="Pearson Correlation Between AIS Variables",
         lower.panel=panel.shade, upper.panel=panel.pie)


####1.d Data Transformations##############################################################
#set Sex as Factor
#Transform SSF to correct for skewedness
ais$Sex=as.factor(ais$Sex)
ais$log.SSF=log(ais$SSF)


#Box Cox transformation analysis for Bfat lambda = 0.711111 (relationship with SSF is strong)
#without SSF in data set lambda = 1
#Opted not to transform Bfat as variables are yet to be determined
library(MASS)
bc = boxcox(Bfat ~.-Sport-SSF, data = ais, lambda = seq(-.2, 2, length = 20))
lambda = bc$x[which.max(bc$y)]
lambda


#######STEP 2 REMOVAL OF VARIABLES#########################################
##########################################################################
#Option 1 Removal of Sport given weak relationship with bodyfat, overall
##Remove Sport from dataset.  No distinct groups when reviewed against bodyfat.  
#Small subsets (gym = 5), uneven variance.  adds risk to fit of model.
lmfit.full=lm(Bfat~.-Sport-SSF, data=ais)
lmfit

#Option 2 removal of variables using stepwise regression
stepfit.full=lm(Bfat~.-Sport-SSF, data=ais)
sfit = step(stepfit.full, direction = "both")
sfit

#Option 3 using regsubsets for variable selection given above assessment show
#Best Subset Selection
library(leaps)
regfit.full=regsubsets(Bfat~.-Sport-SSF, data=ais, method="seqrep", nvmax=11)

plot(regfit.full)
regfit.summary=summary(regfit.full)
plot(regfit.summary$bic, xlab="Number of variables", ylab="BIC", type="l", lwd=2)
minbic=which.min(regfit.summary$bic);minbic
mincp=which.min(regfit.summary$cp);mincp
minadjr2=which.max(regfit.summary$adjr2);minadjr2

#### Function for predictions based on regsubsets objects###
predict.regsubsets <- function(object, newdata, id, ...){
  form <- as.formula(object$call[[2]])
  mat <- model.matrix(form, newdata)
  coefi <-  coef(object, id=id)
  xvars <- names(coefi)
  mat[ , xvars ] %*% coefi
  
} #end function

### Do cross validation with regsubsets model
n = dim(ais)[1] #number of rows in dataset
k = 10 #10-fold cross-validation
groups = c(rep(1:k, floor(n/k)), 1:(n-floor(n/k) * k)) #list of group labels
set.seed(8)
cvgroups = sample(groups,n)
group.error = matrix(nr=11, nc=k) #row = number of variables, column = which fold

for (i in 1:k){
  
  groupi = (cvgroups == i)
  #do regsubsets on the training data
  cv.fit = regsubsets(Bfat~.-Sport-SSF,data=ais[!groupi, ], nvmax=11)
  
  #do predictions for each model size
  for(j in 1:11){
    y.pred = predict.regsubsets(cv.fit, newdata=ais[groupi, ], id=j)
    group.error[j, i] = mean((ais$Bfat[groupi] - y.pred)^2)
  } #end iter over model size
} #end iter over folds

#get the true mean squared error
MSE = apply(group.error, 1, mean)
plot(MSE, xlab="Number of Variables", main="Regsubsets: MSE by Variable Selection")
minmse=which.min(MSE);minmse
MSE[minmse]
se = apply(group.error, 1, sd)/sqrt(k)
se[minmse]

which(MSE <=MSE[minmse]+se[minmse])
minmodel=min(which(MSE <=MSE[minmse]+se[minmse]))
coef(regfit.full, minmodel)
coef(regfit.full, minbic) #coefficients between BIC and MSE within 1 SD are the same.
coef(regfit.full,4) #

#now compare multiple linear models - full, minimum BIC, minimum predictors within one standard deviation of the minimum BIC using anova 
lmfit.full=lm(Bfat~.-Sport-SSF, data=ais)
lmfit.sfit=lm(Bfat~Sex+Ht+Wt+LBM+BMI+log.SSF, data=ais) #subset and lowest bic (6) returned same best subset
lmfit.subsetmin=lm(Bfat~Sex+Wt+Ht, data=ais) #min predictors within 1 standard deviation

anova(lmfit.full, lmfit.sfit) #Min predictors is signficantly different P value 0.8091
anova(lmfit.full, lmfit.subsetmin) #interesting not signficant difference P value 2.2e-16

####Step 3.  Test Proposed Models for heteroskedacity######################################
#Breush Pagan test for heteroscedascity full fit log.SFF - 0.009385
require(lmtest)
bptest(Bfat~.-Sport-SSF, data=ais)

#Breush Pagan test for heteroscedascity model stepwise and regsubsets best bic, cp - selection 0.002447
require(lmtest)
bptest(Bfat~Sex+Ht+Wt+LBM+BMI+log.SSF, data=ais)

#Breush Pagan test for heteroscedasticity 0.06461 #removing variables highly correlated
bptest(Bfat~Sex+Wt+BMI+log.SSF, data=ais)



#######Step 4 Alternative Models to Fit Specifics of the Data##############################
#Random Forest
library(randomForest)
set.seed(8)
ais.rf=randomForest(Bfat~.-Sport-SSF, data=ais, mtry= 11, importance=T)
ais.rf
plot(ais.rf)
#variable importance 
importance(ais.rf)
varImpPlot(ais.rf, main="Random Forest Variable Importance")
ais.all$importance
#tune Random Forest to establish optimal mtry.  Shows default p/3 would likely have higher error
#selected 6 as graph tends to decrease at a signicantly lesser rate between 6-11
x=ais[,c(-1,-2,-13)] #remove Sport, Bfat, and SSF
y=ais[,c(2)] #Bfat
tuneRF(x, y)

####Step 5. Final Models Proposed##########################################################
#Option 1 R2=.9887, AdjR2=.988  All variables
lmfit.full=lm(Bfat~.-Sport-SSF, data=ais)
lmfit.full
summary(lmfit.full)

#Option 2 R2=0.9885 adjR2=0.9882  Regsubsets min BIC, CP and stepwise selection AIC
lmfit.subset=lm(Bfat~Sex+Ht+Wt+LBM+BMI+log.SSF, data=ais)
lmfit.subset
summary(lmfit.subset)

#Option 3 R2=0.6988 adjR2=0.6943  Regsubsets min number of predictors within 1 standard deviation
lmfit.subset2=lm(Bfat~Sex+Ht+Wt, data=ais)
lmfit.subset2
summary(lmfit.subset2)
lmfit.subset2$coefficients


#Option 4 R2=.9459 adjR2=0.9448 Minimum predictors pass Breusch Pagan Test
lmfit.subset3=lm(Bfat~Sex+Wt+BMI+log.SSF, data=ais)
lmfit.subset3
summary(lmfit.subset3)

#Option 5 Random Forest 11 Predictors Considered at each branch
ais.all =randomForest(Bfat~.-Sport-SSF, data=ais, mtry= 11, importance=T)


#Option 6 Random Forest 6 Predictors Considered at each branch
ais.sub=randomForest(Bfat~.-Sport-SSF, data=ais, mtry= 6, importance=T)
ais.sub


#########Double Cross Validation For Model Assessment and Selection########################

#open data
ais=read.csv("ais.csv")

#establish dimensions of dataset - 202 rows, 13 variables
dim(ais)
#Set Sex as Factor
ais$Sex=as.factor(ais$Sex)
ais$log.SSF=log(ais$SSF)

#Rearranged data in the dataframe and removed untransformed variables BMI and SSF
ais.dataframe=data.frame(ais)[,c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14) ]

##### model assessment OUTER shell #####
nvalid = 80 #40% of records
xy.out = ais.dataframe
n.out = dim(xy.out)[1]
#define the validation set 
set.seed(8)
validset = sample(1:n.out,nvalid)
 
trainxy.out = xy.out[-validset,]
testxy.out = xy.out[validset,]
###        inputs trainxy.out       ###
###        :    :    :    :    :    ###
###   entire model-fitting process  ###

##############################
##entire model-fitting process##
xy.in = trainxy.out
n.in = dim(xy.in)[1]
n.in
ncv = 10
if ((n.in%%ncv) == 0) {
  groups.in= rep(1:ncv,floor(n.in/ncv))} else {
    groups.in=c(rep(1:ncv,floor(n.in/ncv)),(1:(n.in%%ncv)))
  }
cvgroups.in = sample(groups.in,n.in)
# with model selection 
allpredictedcv10 = matrix(rep(0,n.in*6),ncol=6)
for (i in 1:ncv) {
  newdata.in = xy.in[cvgroups.in==i,]
  lmfit.all = lm(Bfat~., data=xy.in, subset=(cvgroups.in!=i))
  allpredictedcv10[cvgroups.in==i,1] = predict(lmfit.all,newdata.in)
  
  lmfit.subset1 = lm(Bfat~Sex+Ht+Wt+LBM+BMI+log.SSF, data= xy.in, subset=(cvgroups.in!=i))
  allpredictedcv10[cvgroups.in==i,2] = predict(lmfit.subset1,newdata.in)
  
  lmfit.subset2 = lm(Bfat~Sex+Wt+Ht, data= xy.in, subset=(cvgroups.in!=i))
  allpredictedcv10[cvgroups.in==i,3] = predict(lmfit.subset2,newdata.in)
  
  lmfit.subset3 = lm(Bfat~Sex+Wt+BMI+log.SSF, data= xy.in, subset=(cvgroups.in!=i))
  allpredictedcv10[cvgroups.in==i,4] = predict(lmfit.subset3,newdata.in)
  
  rf.all=randomForest(Bfat~., data=xy.in, mtry= 11, importance=T)
  allpredictedcv10[cvgroups.in==i,5] = predict(rf.all, newdata.in, n.trees=1000, type="response")
  
  rf.sub=randomForest(Bfat~., data=xy.in, mtry= 6, importance=T)
  allpredictedcv10[cvgroups.in==i,6] = predict(rf.sub, newdata.in, n.trees=1000, type="response")

}

allmodelcv.in = rep(0,6)
for (m in 1:6) 
{allmodelcv.in[m] = sum((allpredictedcv10[,m]-xy.in$Bfat)^2)/n.in}
plot(1:6,allmodelcv.in,col="red",pch=20, xlab="Model", ylab="Model CV Error (Training Data)")  
legend("topleft", legend=c(paste("Model 1 ",format(allmodelcv.in[1], digits=4)), paste("Model 2 ",format(allmodelcv.in[2], digits=4)),
paste("Model 3 ",prettyNum(allmodelcv.in[3], digits=4)), paste("Model 4 ",prettyNum(allmodelcv.in[4], digits=4)),
paste("Model 5 ",format(allmodelcv.in[5], digits=4)), paste("Model 6 ",format(allmodelcv.in[6], digits=4)))) # paste("Model 2",format(allmodelcv.in[2], digits=4))
bestmodel.in = (1:6)[order(allmodelcv.in)[1]]  # actual selection
###      resulting in bestmodels     ###
bestmodel = ifelse(length(bestmodel.in)==1,bestmodel.in,sample(bestmodel.in,1))

if (bestmodel == 1)  {
  lmfit.all.train = lm(Bfat~., data=trainxy.out)
  predictvalid = predict(lmfit.all.train, testxy.out)
  lmfit.all.full = lm(Bfat~., data=xy.out)
  predictfull=predict(lmfit.all.full, xy.out)
  bestmodelformula=lmfit.all.full
}

if (bestmodel == 2)  {
  lmfit.subset1.train = lm(Bfat~Sex+Ht+Wt+LBM+BMI+log.SSF, data=trainxy.out)
  predictvalid = predict(lmfit.subset1.train, testxy.out)
  lmfit.subset1.full = lm(Bfat~Sex+Ht+Wt+LBM+BMI+log.SSF, data=xy.out)
  predictfull = predict(lmfit.subset1.full, xy.out) 
  bestmodelformula=lmfit.subset1.full
}
if (bestmodel == 3)  {
  lmfit.subset2.train = lm(Bfat~Sex+Wt+Ht, data=trainxy.out)
  predictvalid = predict(lmfit.subset2.train, testxy.out)
  lmfit.subset2.full = lm(Bfat~Sex+Wt+Ht, data=xy.out)
  predictvalid = predict(lmfit.subset2.full, xy.out)
  bestmodelformula=lmfit.subset2.full
}
if (bestmodel == 4)  {
  lmfit.subset3.train = lm(Bfat~Sex+Wt+BMI+log.SSF, data=trainxy.out)
  predictvalid = predict(lmfit.subset3.train, testxy.out)
  lmfit.subset3.full = lm(Bfat~Sex+Wt+BMI+log.SSF, data=xy.out)
  predictfull = predict(lmfit.subset3.full, xy.out)
  bestmodelformula=lmfit.subset3.full
}
if (bestmodel == 5)  {
  rf.all.train = randomForest(Bfat ~ ., data=trainxy.out, mtry= 11, importance=T )
  predictvalid = predict(rf.all.train, testxy.out, n.trees=1000, type="response")
  rf.all.full = randomForest(Bfat ~ ., data=xy.out, mtry= 11, importance=T )
  predictfull = predict(rf.all.full, xy.out, n.trees=1000, type="response")
  bestmodelformula=rf.all.full
  
}

if (bestmodel == 6)  {
   rf.sub.train = randomForest(Bfat ~ ., data=trainxy.out, mtry=6 , importance=T )
   predictvalid = predict(rf.sub, testxy.out, n.trees=1000, type="response")
   rf.sub.full = randomForest(Bfat ~ ., data=xy.out, mtry=6 , importance=T )
   predictfull = predict(rf.sub, xy.out, n.trees=1000, type="response")
   bestmodelformula=rf.sub.full
}
allmodelcv.in
#assessment using validation set
testy.out = testxy.out$Bfat ### or is this testxy.out?
validCV.out = sum((predictvalid-testy.out)^2)/nvalid; validCV.out
validR2.out = 1-sum((predictvalid-testy.out)^2)/sum((testy.out-mean(testy.out))^2); validR2.out


######Best Model Fit on Entire Data Set###############
bestmodel
y.out=xy.out$Bfat
CV.out = sum((predictfull-y.out)^2)/n.out; CV.out
R2.out = 1-sum((predictfull-y.out)^2)/sum((y.out-mean(y.out))^2); validR2.out
plot(bestmodelformula, main="Best Model Error")




 
















