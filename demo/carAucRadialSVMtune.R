########## R script: carAucRadialSVMtune ##########

# For doing tuning a radial basis function support vector 
# machines for classification based on car auction data
# from the Kaggle platform (www.kaggle.com) ``Don't Get 
# Kicked!'' competition (held during 2012).

# Last changed: 20 JUN 2017

# Set main variables and seed:

indctrInds <- c(3,6,10,13,29,30,32,37,40,41,43)
contnsInds <- c(9,31,39,47)
set.seed(1)

# Load required packages:

library(e1071) ; library(HRW)

# Read in original training data:

data(carAuction)

# Create new training data such that the
# number of "good buy" and "bad buy" cars
# are equal.

nOrig <- nrow(carAuction)
indsBadBuy <- (1:nOrig)[carAuction[,2] == 1]
indsGoodBuy <- setdiff(1:nOrig,indsBadBuy)
indsGoodBuySub <- sample(indsGoodBuy,length(indsBadBuy),replace = FALSE)
allData <- carAuction[c(indsGoodBuySub,indsBadBuy),]

# Break off a random subset of size 1000 to use for estimating
# test error, with remaining data used for training.

testInds <- sample(1:nrow(allData),1000,replace = FALSE)
trainInds <- setdiff(1:nrow(allData),testInds)
trainData <- allData[trainInds,]

# Obtain a subset of size 1000 for tuning
# (otherwise the tuning can be very slow):

indsTune <- sample(1:nrow(trainData),1000,replace = FALSE)
trainData <- trainData[indsTune,]

# Set response variable (y = 1 if bad buy and y = -1 if good buy) 
# and the the matrix containing all possibly predictors

yTrain <- 2*trainData[,2] - 1
predsAllTrain <-  trainData[,-c(1,2)]
predMatTrain <- cbind(predsAllTrain[,indctrInds],predsAllTrain[,contnsInds])
trainDF <- data.frame(yTrain,predMatTrain)

tune.svmObj <- tune.svm(factor(yTrain)~.,data = trainDF,gamma = 10^(-5:-1),
                        cost = 10^seq(1,2,length = 5))
      
cat("The selected radial basis function SVM parameters are:\n")
print(tune.svmObj$best.parameters)

############ End of carAucRadialSVMtune ############

