########## R script: carAucPenSplSVMtune ##########

# For tuning penalized spline support vector machines
# for classification based on the car auction data.

# Last changed: 03 AUG 2017

# Set main variables and seed:

indctrInds <- c(3,6,10,13,29,30,32,37,40,41,43)
contnsInds <- c(9,31,39,47)
set.seed(1)

# Set cross-validation size and grid of lambda values:

kCVval <- 2
lambdaLow <- 10 ; lambdaUpp <- 100
numLambda <- 21
lambdaGrid <- exp(seq(log(lambdaLow),log(lambdaUpp),length = numLambda))

# Load required packages:

library(LowRankQP)
library(HRW)

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

y <- 2*trainData[,2] - 1
predsAll <-  trainData[,-c(1,2)]
n <- length(y)  
numLin <- length(indctrInds)   ;   numNonLin <- length(contnsInds)

# Obtain vector of names of chosen predictors:

namesVec <- names(predsAll)[c(indctrInds,contnsInds)]

# Obtain X matrix in original units:

Xorig <- cbind(rep(1,n),predsAll[,indctrInds],predsAll[,contnsInds])

# Transform all non-intercept columns to the unit interval
# for penalized spline support vector machine fitting:

origMins <- apply(Xorig,2,min) ; origMaxs <- apply(Xorig,2,max) 
X <- matrix(NA,nrow(Xorig),ncol(Xorig))
X[,1] <- 1
for (j in 2:ncol(Xorig))
   X[,j] <- (Xorig[,j] - origMins[j])/(origMaxs[j] - origMins[j])
ncX <- ncol(X)

nVal <- nrow(trainData)
scramInds <- sample(1:nVal,nVal,replace = FALSE)
bdyInds <- round(seq(1,nVal,length = (kCVval+1)))
sttInds <- bdyInds[-(kCVval+1)]
endInds <- c(sttInds[2:kCVval] + 1,nVal)
CVinds <- vector("list",kCVval)
for (iCV in 1:kCVval)
   CVinds[[iCV]] <- scramInds[sttInds[iCV]:endInds[iCV]]

misClassMat <- matrix(NA,kCVval,numLambda)

for (iCV in 1:kCVval)
{
   testInds <- CVinds[[iCV]]
   trainInds <- setdiff(1:nVal,testInds)
   yTest <- y[testInds]
   yTrain <- y[trainInds]
   nTrain <- length(yTrain)

   Xtest <- X[testInds,]
   Xtrain <- X[trainInds,]
   
   intKnots <- vector("list",numNonLin) 
   numIntKnots <- c(5,rep(15,3))
   Ztrain <- NULL ; Ztest <- NULL
   for (ell in 1:numNonLin)
   {
      j <- ell + numLin
      xTrainCurr <- Xtrain[,j]
      xTestCurr <- Xtest[,j]
      intKnotsCurr <- quantile(unique(xTrainCurr),seq(0,1,length = numIntKnots[ell]+2)
                              [-c(1,numIntKnots[ell]+2)])
      Ztrain <- cbind(Ztrain,ZOSull(xTrainCurr,c(0,1),intKnotsCurr))
      Ztest <- cbind(Ztest,ZOSull(xTestCurr,c(0,1),intKnotsCurr))
      intKnots[[ell]] <- intKnotsCurr
   }
      
   yXtrain <- as.matrix(yTrain*Xtrain) 

   for (iLambda in 1:numLambda)
   {
      lambda <- lambdaGrid[iLambda]

      ZtrainTilde <- Ztrain*(1/sqrt(2*lambda))
      yZtrainTilde <- as.matrix(yTrain*ZtrainTilde)

      LowRankQPobj <- LowRankQP(yZtrainTilde,rep(-1,nTrain),
                                t(yXtrain),rep(0,ncol(Xtrain)),rep(1,nTrain))

      # Obtain estimated beta and u vectors:

      betaHat <- as.vector(LowRankQPobj$beta) 
      uHat <-  as.vector(crossprod(yZtrainTilde,LowRankQPobj$alpha))

      ZtestTilde <- Ztest*(1/sqrt(2*lambda))
      Xtest <- as.matrix(Xtest)

      fHatSVM <- Xtest%*%betaHat + ZtestTilde%*%uHat
      classVec <- sign(fHatSVM)
      confMat <- as.matrix(table(yTest,classVec))
 
      misClassMat[iCV,iLambda] <- (100*(sum(confMat) 
                                         - sum(diag(confMat)))/sum(confMat))
   }
}

misClassMed <- apply(misClassMat,2,median)
plot(lambdaGrid,misClassMed,type = "n",bty = "l",xlab = expression(lambda),
     ylab = "average 2-fold cross-validation misclassification rate")
lines(lambdaGrid,misClassMed,col = "blue",lwd = 2)

############ End of carAucPenSplSVMtune ############

