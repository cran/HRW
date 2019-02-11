########## R script: carAucRadialSVM ##########

# For fitting and evaluating a radial basis function
# support vector machines for classification based on
# the car auction data from the Kaggle platform
# (www.kaggle.com) ``Don't Get  Kicked!'' competition
# (held during 2012).

# Last changed: 20 JUN 2017

# Set main variables and seed:

indctrInds <- c(3,6,10,13,29,30,32,37,40,41,43)
contnsInds <- c(9,31,39,47)
set.seed(1)

# Load required packages:

library(e1071) ; library(HRW)

# Load in original training data:

data(carAuction)

# Create new training data such that the
# numbers of "good buy" and "bad buy" cars
# are equal.

nOrig <- nrow(carAuction)
indsBadBuy <- (1:nOrig)[carAuction[,2] == 1]
indsGoodBuy <- setdiff(1:nOrig,indsBadBuy)
indsGoodBuySub <- sample(indsGoodBuy,length(indsBadBuy),
                          replace = FALSE)
allData <- carAuction[c(indsGoodBuySub,indsBadBuy),]

# Break off a random subset of size 1000 to use for estimating
# test error, with remaining data used for training. 

testInds <- sample(1:nrow(allData),1000,replace = FALSE)
trainInds <- setdiff(1:nrow(allData),testInds)
testData <-  allData[testInds,]
trainData <- allData[trainInds,]

# Set response variable (y = 1 if bad buy and y = -1 if good buy) 
# and the the matrix containing all possibly predictors

yTrain <- 2*trainData[,2] - 1
predsAllTrain <-  trainData[,-c(1,2)]
predMatTrain <- cbind(predsAllTrain[,indctrInds],predsAllTrain[,contnsInds])
trainDF <- data.frame(yTrain,predMatTrain)

# Set gamma and cost values to be those from tuning:

gammaVal <- 0.01 ; costVal <- 17.78279
fitSVM <- svm(factor(yTrain)~.,data = trainDF,gamma = gammaVal,cost = costVal)   

yTest <- testData[,2]
predMatTest <-  testData[,-c(1,2)]

# Obtain classifications based on test data:

yHat <- 2*as.numeric(predict(fitSVM,predMatTest)) - 3

# Obtain and print out the confusion matrix and test error:

confMat <- as.matrix(table(yTest,yHat))

names(dimnames(confMat)) <- rep("",2)
dimnames(confMat)[[1]] <- c("really good buy","really bad buy")
dimnames(confMat)[[2]] <- c("classified good buy","classified bad buy")

cat("\n")
cat("The confusion matrix is:\n")
print(confMat)

cat("\n")
misClassError <- round(100*(sum(confMat) - sum(diag(confMat)))/sum(confMat),3)

# Obtain and estimate misclassifcation rate:

misClassRate <- (100*(sum(confMat) - sum(diag(confMat)))/sum(confMat))

cat("The estimated misclassification rate is: ",misClassRate,"%\n\n",sep = "")

############ End of carAucRadialSVM ############

