#############Step 1 Read in Data and Analyze###############################
library(stats)
library(pROC)
library(tree)
#Read in CTG Raw Data
ctg_raw=read.csv('CTG_RAW.csv')

#Get dimensions of ctg_raw dataframe
dim(ctg_raw)
#Review column names
colnames(ctg_raw)
#review first few rows
head(ctg_raw)

#summarize dataset - look for NA's and columns with no values
summary(ctg_raw)

#remove blank rows and summary counts at end of file (lines 2, 2128-2130)
ctg_raw=na.omit(ctg_raw)

#re-review dataframe
summary(ctg_raw)


#create data frame with results for comparison.  
#Removed descriptor columns: FileName, Date, b (start instant), e(end instant), LBE (baseline Medical expert).  Kept LB (sisPorto)
#Removed classification columns used for secondary classification "CLASS" as these are an alternative classification method not explored here.
#A, B,C, D, E, AD, DE, LD, FS, SUSP, CLASS
#Removed columns with no values DR  (repetitive decelerations = 0 for all rows). 

#Dataframe with features and one response variable NSP
ctg_raw_baseline=ctg_raw[ ,c(7:17, 19:28, 40)]

#Perform column sums to ensure all numeric values
colSums(ctg_raw_baseline)

###############Histograms###########################
###is it a gaussian distribution?
attach(ctg_raw_baseline)

hist(LB)      #normal
hist(AC)      #right skewed with one dominant category.  long tail
hist(FM)      #right skewed with one dominant category.  long tail
hist(UC)      #right skewed
hist(ASTV)    #potentially two peaks
hist(MSTV)    #right skewed 
hist(ALTV)    #right skewed - one dominant category with long tail
hist(MLTV)    #right skewed
hist(DL)      #rightskewed Poisson (0 is dominant category with long tail (1-12))
hist(DS)      #right skewed Poisson (0, with 1 being an outlier)
hist(DP)      #right skewed Poisson (0, 1, 2,3)
hist(Width)   #potentially two peaks
hist(Min)     #potentially two peaks
hist(Max)     #right skewed.  generally  normal.  long right tail
hist(Nmax)    #left skewed Poisson (0~20)
hist(Nzeros)  # lef skewed Poisson (0,1,2,3)  
hist(Mode)    #slightly left skewed. generally normal
hist(Mean)    #slightly left skewed. generally normal
hist(Median)   #normal
hist(Variance) #right skewed
hist(Tendency) #Poisson (-1, 0, 1)


######Perform PCA Analysis on Each Response Variable to See Distribution################
###################create data frame for predictions without results.  Features only.
ctg_raw_predictions=ctg_raw_baseline[,c(-22)] #remove response variable
#look for missing values
summary(ctg_raw_predictions)

#perform log transformations based on histograms
ctg_raw_predictions_log=log(ctg_raw_predictions+1.01) #need to add 1.01 to deal with -1 values of tendency and other variable producing INF
 

#look for NA's
summary(ctg_raw_predictions_log)
# include scaling due to widely varied magnitudes
ctg.pc.info = prcomp(ctg_raw_predictions,center=T,scale=T)
#review comparison to log transformation
ctg.pc.info2 = prcomp(ctg_raw_predictions_log,center=T,scale=T)

ctg.pc.info$rotation[,1:2]  #loadings
ctg.pc.info2$rotation[,1:2]
###no log transformations needed.  no multiplicatively additive correlation

par(mfrow=c(1,1))
my.col=rep("black", 2126)
my.col[which(ctg_raw_baseline$NSP=="1")]="blue"
my.col[which(ctg_raw_baseline$NSP=="2")]="red"
my.pch=rep(17, 2126)
my.pch[which(ctg_raw_baseline$NSP=="1")] = 1
my.pch[which(ctg_raw_baseline$NSP=="2")] = 19
plot(ctg.pc.info$x[,1:2],col=my.col, pch=my.pch, main="Physician Assessed NSP Levels")
legend("topright", legend=c("Normal", "Suspect", "Pathologic"), col=c("blue", "red", "black"),  pch=c(1, 19, 17))
plot(ctg.pc.info2$x[,1:2],col=my.col, pch=my.pch, main="log")
legend("topright", legend=c("Normal", "Suspect", "Pathologic"), col=c("blue", "red", "black"),  pch=c(1, 19, 17))

#check number of PC's to use
par(mfrow=c(1,1))
plot(ctg.pc.info, main="Importance of Components All Responses")
summary(ctg.pc.info)

#calculate PVE .95 at PC14, .99 PC17 (PC20, PC21 at 1)
ctg.pc.info$sdev
vjs = ctg.pc.info$sdev^2
pve = vjs/sum(vjs)
cumsum(pve)

#Plot of cumulative sums - scree plot
plot(cumsum(pve), type = "o", ylab="Cumulative PVE", xlab="Principal Component")

#which variables contribute the most.  PC1 and PC2 contribute the most
#arrows pointing in opposite directions and of differing lengths in combined set
#
#all variables are contributors it seems, but of oppposite directions and weights
par(mfrow=c(1,1))
biplot(ctg.pc.info,scale=0, main="All Results")

ctg.pc.info$rotation[,1]  # loadings for first principal component
ctg.pc.info$rotation[,2]  # loadings for second principal component
ctg.pc1scores = ctg.pc.info$x[,1]  # first principal component score vector 
#pc1scores
ctg.pc2scores = ctg.pc.info$x[,2]  # second principal component score vector


####################Analyze just the NSP level 1 (Normal)###################
ctg_raw_limited_1=ctg_raw_baseline[ctg_raw_baseline$NSP=="1", ] ###level 
ctg_raw_limited_predictions_1=ctg_raw_limited_1[,c(-22)] #has DS, Removed NSP
#Review for values in DS Feature
summary(ctg_raw_limited_predictions_1)

# include scaling and centering due to widely varied magnitudes and distribution
ctg.ltd1.pc.info = prcomp(ctg_raw_limited_predictions_1,center=T,scale=T)
ctg.ltd1.pc.info$rotation  #loadings

#check number of PC's to use
par(mfrow=c(1,1))
plot(ctg.ltd1.pc.info, main="Importance of Components For Normal (NSP=1)")
summary(ctg.ltd1.pc.info)


#calculate PVE .95 at PC14, .99 PC17 (PC20, PC21 at 1)
ctg.ltd1.pc.info$sdev
vjs1 = ctg.ltd1.pc.info$sdev^2
pve1 = vjs1/sum(vjs1)
cumsum(pve1)

#Plot of cumulative sums - scree plot
plot(cumsum(pve1), type = "o", ylab="Cumulative PVE", xlab="Principal Component", main="Normal (NSP=1)")

#which variables contribute the most.  PC1 and PC2 distinctly contribute most
#arrows in biplot indicate a strong movement in the same direction
#some features stronger than others based on arrow lenght
#PC20, 21 may not be needed given the variance is explained 99% at PC19
par(mfrow=c(1,1)) 
biplot(ctg.ltd1.pc.info,scale=0, main="Normal (NSP=1)")
ctg.ltd1.pc.info$rotation[,1]  # loadings for first principal component
ctg.ltd1.pc.info$rotation[,2]  # loadings for second principal component
ctg.ltd1.pc1scores = ctg.ltd1.pc.info$x[,1]  # first principal component score vector 
ctg.ltd1.pc2scores = ctg.ltd1.pc.info$x[,2]  # second principal component score vector

####################Analyze just the NSP level 3 (Pathologic)###################

ctg_raw_limited_3=ctg_raw_baseline[ctg_raw_baseline$NSP=="3", ] ###level 
ctg_raw_limited_predictions_3=ctg_raw_limited_3[,c(-22)] #excludes NSP
#Review for category with missing values
summary(ctg_raw_limited_predictions_3)
# include scaling and centering due to widely varied magnitudes and non-normal distributions
ctg.ltd3.pc.info = prcomp(ctg_raw_limited_predictions_3,center=T,scale=T)
ctg.ltd3.pc.info$rotation  #loadings

#check number of PC's to use
par(mfrow=c(1,1))
plot(ctg.ltd3.pc.info, main="Importance of Components For Pathologic (NSP=3)")
summary(ctg.ltd3.pc.info)


#calculate PVE .95 at PC12, .99 PC17 (PC20, PC21 at 1)
ctg.ltd3.pc.info$sdev
vjs3 = ctg.ltd3.pc.info$sdev^2
pve3 = vjs3/sum(vjs3)
cumsum(pve3)

#Plot of cumulative sums - scree plot
plot(cumsum(pve3), type = "o", ylab="Cumulative PVE", xlab="Principal Component", main="Pathologic (NSP=3)")

#which variables contribute the most.  only show the first two (highest)
#look at arrows and length (same direction, 
#Two distinct groups in distinctly oppositive directions
#all variables are contributors it seems.  Variance explained 99% at PC17.  May not need 18-21PC components 20-21 could be excluded
par(mfrow=c(1,1))
biplot(ctg.ltd3.pc.info,scale=0, main="Pathologic (NSP=3)")

ctg.ltd3.pc.info$rotation[,1]  # loadings for first principal component
ctg.ltd3.pc.info$rotation[,2]  # loadings for second principal component
ctg.ltd3.pc1scores = ctg.ltd3.pc.info$x[,1]  # first principal component score vector 
ctg.ltd3.pc2scores = ctg.ltd3.pc.info$x[,2]  # second principal component score vector

####################Analyze just the NSP level 2 (Suspect)###################

ctg_raw_limited_2=ctg_raw_baseline[ctg_raw_baseline$NSP=="2", ] ###level
ctg_raw_limited_predictions_2=ctg_raw_limited_2[,c(-10,-22)]  #removed DS all zeros

summary(ctg_raw_limited_predictions_2)


# include scaling due to widely varied magnitudes
ctg.ltd2.pc.info = prcomp(ctg_raw_limited_predictions_2,center=T,scale=T)
ctg.ltd2.pc.info$rotation  #loadings

#check number of PC's to use
par(mfrow=c(1,1))
plot(ctg.ltd2.pc.info, main="Importance of Components For Pathologic (NSP=3)")
summary(ctg.ltd2.pc.info)

#calculate PVE .95 at PC12, .99 PC16 (PC19, PC20 at 1) - due to reduction of DS
ctg.ltd2.pc.info$sdev
vjs2 = ctg.ltd2.pc.info$sdev^2
pve2 = vjs2/sum(vjs2)
cumsum(pve2)

#Plot of cumulative sums - scree plot
plot(cumsum(pve2), type = "o", ylab="Cumulative PVE", xlab="Principal Component", main="Suspect (NSP=2)")

#PC1 and PC2 highest.  PC3 also adds
#tightly grouped with left leaning outliers
#majority of arrows leaning in similiar direction to normal plot exceptions are Min, ALTV, ASTV
#all components contribute.  may drop PC's 20-21
biplot(ctg.ltd2.pc.info,scale=0, main="Suspect without DS (NSP=2) ")
ctg.ltd2.pc.info$rotation[,1]  # loadings for first principal component
ctg.ltd2.pc.info$rotation[,2]  # loadings for second principal component
ctg.ltd2.pc1scores = ctg.ltd2.pc.info$x[,1]  # first principal component score vector 
ctg.ltd2.pc2scores = ctg.ltd2.pc.info$x[,2]  # second principal component score vector

############compare 3 biplots###################
par(mfrow=c(1,3)) 
biplot(ctg.ltd1.pc.info,scale=0, main="Normal (NSP=1)")
biplot(ctg.ltd2.pc.info,scale=0, main="Suspect without DS (NSP=2)")
biplot(ctg.ltd3.pc.info,scale=0, main="Pathologic (NSP=3)")
par(mfrow=c(1,1))
biplot(ctg.pc.info,scale=0, main="All Results")

ctg.ltd1.pc.info$rotation[,1]  # loadings for first principal component
ctg.ltd1.pc.info$rotation[,2]
ctg.ltd3.pc.info$rotation[,1]  # loadings for first principal component
ctg.ltd3.pc.info$rotation[,2]


###################Analyze just the NSP level 1 and 3 (Normal and Pathologic) ###################

ctg_raw_limited_13=ctg_raw_baseline[ctg_raw_baseline$NSP!="2", ] ###level 
ctg_raw_limited_predictions_13=ctg_raw_limited_13[,c(-22)] #has DS
 
#set response as factor for predictions later
ctg_raw_limited_13$NSP=as.factor(as.numeric(ctg_raw_limited_13$NSP))
#set levels to match NSP levels for confusion matrix later
levels(ctg_raw_limited_13$NSP)= c(1, 3)
#validate dataframe
summary(ctg_raw_limited_13$NSP)
summary(ctg_raw_limited_predictions_13)

# include scaling and centering due to widely varied magnitudes & distributions
ctg.ltd13.pc.info = prcomp(ctg_raw_limited_predictions_13,center=T,scale=T)
ctg.ltd13.pc.info$rotation  #loadings
ctg.ltd13.pc.info$x

#check number of PC's to use
par(mfrow=c(1,1))
plot(ctg.ltd13.pc.info, main="Importance of Components For Normal (NSP=1) and Pathologic (NSP=3)")
summary(ctg.ltd13.pc.info)


#calculate PVE .95 at PC14, .99 PC17 (PC20, PC21 at 1)
ctg.ltd13.pc.info$sdev
vjs13 = ctg.ltd13.pc.info$sdev^2
pve13 = vjs13/sum(vjs13)
cumsum(pve13)

#Plot of cumulative sums - scree plot
plot(cumsum(pve13), type = "o", ylab="Cumulative PVE", xlab="Principal Component", main="NSP 1 and 3")
plot(pve13, type = "b", ylab="Proporation Variance Explained", xlab="Principal Component", main="NSP 1 and 3")

#which variables contribute the most.  only show the first two (highest)
#look at arrows and length (same direction, contributing similiar amounts in the frist two PC)
#loadings tell us how much they move in each direction
#all variables are contributors it seems
par(mfrow=c(1,1)) 
biplot(ctg.ltd13.pc.info,scale=0, main="Normal and Pathologic (NSP=1 or NSP=3)")

ctg.ltd13.pc.info$rotation[,1]  # loadings for first principal component
ctg.ltd13.pc.info$rotation[,2]  # loadings for second principal component
ctg.ltd1.pc1scores = ctg.ltd1.pc.info$x[,1]  # first principal component score vector 
#pc1scores
ctg.ltd1.pc2scores = ctg.ltd1.pc.info$x[,2]  # second principal component score vector

################Use Principal Components from 1 and 3 in supervised decision tree model################
##############Train and Test on 13#######################################
#Add training set which is all records (NSP=1 or 3) with principal components
train.data=data.frame(NSP=ctg_raw_limited_13$NSP, ctg.ltd13.pc.info$x)

##################Model Assessment - Optimal Leaves##############################

#decision tree model off of training data to determine optimal leaves
set.seed(12)
tree.model=tree(NSP~., data=train.data)
par(mfrow=c(1,1))
plot(tree.model, main="Full Training (NSP=1 and NSP=3)")

tree.predictions=predict(tree.model, train.data, type="class")

##confusion matrix
tree.predict.table=table(tree.predictions,ctg_raw_limited_13$NSP)
tree.predict.table
##Error Rate
(tree.predict.table[1,1]+tree.predict.table[2,2])/(tree.predict.table[1,1]+tree.predict.table[2,2]+tree.predict.table[1,2]+tree.predict.table[2,1])
#true positive rate train2
tree.predict.table[2,2]/(tree.predict.table[2,2] + tree.predict.table[1,2])
#false positive (1-specificity) train2
tree.predict.table[2,1]/(tree.predict.table[2,1] + tree.predict.table[1,1])

#auc
tree.auc=roc(response=tree.predictions, predictor=as.numeric(ctg_raw_limited_13$NSP))
tree.auc$auc


###choose better model through pruning

tree.model.cv=cv.tree(tree.model, K=10, FUN=prune.misclass)
 
plot(tree.model.cv)
min(tree.model.cv$dev)
which(tree.model.cv$dev==min(tree.model.cv$dev))
max(tree.model.cv$size) #max model size
tree.model.cv$size[which(tree.model.cv$dev==min(tree.model.cv$dev))] #15

prune.tree.model15=prune.misclass(tree.model, best=15)
prune.tree.model18=prune.misclass(tree.model, best=18)
 
par(mfrow=c(1,2))
plot(prune.tree.model18, main="Max Model Size 17")
text(prune.tree.model18, pretty=0)
plot(prune.tree.model15, main="Best Model 15")
text(prune.tree.model15, pretty=0)

####predictions on pruned tree###################
prune.tree.predictions=predict(prune.tree.model15, train.data, type="class")

##confusion matrix
prune.tree.predict.table=table(prune.tree.predictions,ctg_raw_limited_13$NSP)
prune.tree.predict.table
##Error Rate (Success)
(prune.tree.predict.table[1,1]+prune.tree.predict.table[2,2])/(prune.tree.predict.table[1,1]+prune.tree.predict.table[2,2]+prune.tree.predict.table[1,2]+prune.tree.predict.table[2,1])
#true positive rate train2
prune.tree.predict.table[2,2]/(prune.tree.predict.table[2,2] + prune.tree.predict.table[1,2])
#false positive (1-specificity) train2
prune.tree.predict.table[2,1]/(prune.tree.predict.table[2,1] + prune.tree.predict.table[1,1])

#auc
prune.auc=roc(response=prune.tree.predictions, predictor=as.numeric(ctg_raw_limited_13$NSP))
prune.auc$auc

####################Cross Validation for Model Assessment###################################
n=dim(train.data)[1]
k=10
set.seed(12)
groups= c(rep(1:k, floor(n/k)), 1:(n-floor(n/k)*k))
cvgroups=sample(groups,n)
CV10predictions=rep(NA,n)
CV10errorrate=rep(NA,k)
CV10truepositive=rep(NA,k)
CV10falsepositive=rep(NA,k)
CV10auc=rep(NA,k)

for (i in 1:k) {
  groupi=(cvgroups==i)  
  train.pca=train.data[-groupi,]
  valid.pca=train.data[groupi, ]
  valid.pca.y=train.data[groupi, 1]
  tree.model=tree(NSP~., data=train.pca)
  tree.predictions.cv=predict(prune.tree.model15, valid.pca, type="class")
  ##confusion matrix
  tree.predict.table=table(tree.predictions.cv,valid.pca.y)
  CV10auc[i]=roc(response=tree.predictions.cv, predictor=as.numeric(valid.pca.y))$auc
  CV10predictions[groupi]=tree.predictions.cv
  CV10errorrate[i]=(tree.predict.table[1,1]+tree.predict.table[2,2])/(tree.predict.table[1,1]+tree.predict.table[2,2]+tree.predict.table[1,2]+tree.predict.table[2,1])
  CV10truepositive[i]=tree.predict.table[2,2]/(tree.predict.table[2,2] + tree.predict.table[1,2])
  CV10falsepositive[i]=tree.predict.table[2,1]/(tree.predict.table[2,1] + tree.predict.table[1,1])
}
which(is.na(CV10predictions)) #ensure all predictions populated (no NA)
mean(CV10truepositive)
mean(CV10errorrate) #success 
mean(CV10falsepositive)
mean(CV10auc)


##################Classify NSP 2 Data##########################
#transform new data into PCA
data.2=ctg_raw_limited_2[,-22]


#PCA predictions to get PC loadings
new.data = predict(ctg.ltd13.pc.info, newdata = data.2)
new.data = as.data.frame(new.data)
new.data


#predict classification using pruned tree model 15
new.data.predictions=predict(prune.tree.model15, newdata=new.data, type="class")
table(new.data.predictions) 


######scatterplot all predictions ###
#transform new data into PCA### 
all.pred.data=predict(ctg.ltd13.pc.info, newdata=ctg_raw_predictions)
all.pred.data=as.data.frame(all.pred.data)
all.pred.data.y=ctg_raw_baseline

all.prune.tree.predictions=predict(prune.tree.model15, newdata=all.pred.data, type="class")
table(all.prune.tree.predictions)
summary(prune.tree.model15)  ###understates as we don't know true outcomes of suspect readings


###plot predictions###############################
#png(file="Model_MD_Compare.png") #image in executive summary
par(mfrow=c(1,2))
my.col=rep("black", 2126)
my.col[which(all.prune.tree.predictions=="1")]="blue"
my.pch=rep(17, 2126)
my.pch[which(all.prune.tree.predictions=="1")] = 1
plot(ctg.pc.info$x[,1:2],col=my.col, pch=my.pch, main="Model Predictions")
legend("topright", legend=c("Normal", "Pathologic"), col=c("blue", "black"),  pch=c(1, 17))

###original plot of physician NSP levels##########################
my.col=rep("black", 2126)
my.col[which(ctg_raw_baseline$NSP=="1")]="blue"
my.col[which(ctg_raw_baseline$NSP=="2")]="red"
my.pch=rep(17, 2126)
my.pch[which(ctg_raw_baseline$NSP=="1")] = 1
my.pch[which(ctg_raw_baseline$NSP=="2")] = 19
plot(ctg.pc.info$x[,1:2],col=my.col, pch=my.pch, main="Physician Assessed NSP Levels")
legend("topright", legend=c("Normal", "Suspect", "Pathologic"), col=c("blue", "red", "black"),  pch=c(1, 19, 17))
#dev.off()



