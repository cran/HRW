########## R script: carAucPenSplSVM ##########

# For fitting and evaluating a penalized 
# spline support vector machines for classification 
# based on the car auction from the Kaggle platform
# (www.kaggle.com) ``Don't Get  Kicked!'' competition
# (held during 2011-2012).

# Last changed: 19 JUN 2017

# Set main variables and seed:

indctrInds <- c(3,6,10,13,29,30,32,37,40,41,43)
contnsInds <- c(9,31,39,47)
set.seed(1)
lambda <- 44.7

# Load required packages:

library(LowRankQP) ; library(HRW)

# Load in original training data:

data(carAuction)

# Set long names vector:

longPredNames <- c("purchased in 2010","auctioned by Adesa","auctioned by Manheim",
                   "2003 vehicle","2004 vehicle","2005 vehicle","2006 vehicle",
                   "2007 vehicle","age at sale","Chevrolet","Ford","Dodge","Chrysler",
                   "trim is Bas","trim is LS","trim is SE","sub-model is 4DSEDAN",
                   "sub-model is 4DESANSL","sub-model is 4DSESANSE",
                   "coloured silver","coloured white","coloured blue",
                   "coloured grey","coloured black","coloured red",
                   "coloured gold","coloured orange","manual transmission",
                   "alloy wheels","covered wheels","odometer reading",
                   "made in U.S.A.","made in Asia (not Japan or South Korea)",
                   "truck","medium-sized","sports utility vehicle","compact","van",
                   "acquisition price","purchased in Texas","purchased in Florida",
                   "purchased in California","purchased in North Carolina",
                   "purchased in Arizona","purchased in Colorado",
                   "purchased in South Carolina","acquisition cost",
                   "online sale","warranty cost")

# Create new training data such that the numbers
# of "good buy" cars and "bad buy" cars are equal.

nOrig <- nrow(carAuction)
indsBadBuy <- (1:nOrig)[carAuction[,2] == 1]
indsGoodBuy <- setdiff(1:nOrig,indsBadBuy)
indsGoodBuySub <- sample(indsGoodBuy,length(indsBadBuy),replace = FALSE)
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
numIndctr <- length(indctrInds)   
numContns <- length(contnsInds)

# Obtain vector of names of chosen predictors:

namesVec <- longPredNames[c(indctrInds,contnsInds)]

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

intKnots <- vector("list",numContns) 
reBlockInds <-  vector("list",numContns) 
currReBlockStt <- 1
numIntKnots <- c(5,rep(15,3))
Z <- NULL
for (ell in 1:numContns)
{
   j <- ell + numIndctr
   xCurr <- X[,j]
   intKnotsCurr <- quantile(unique(xCurr),
                            seq(0,1,length = numIntKnots[ell]+2)
                            [-c(1,numIntKnots[ell]+2)])
   Z <- cbind(Z,ZOSull(xCurr,c(0,1),intKnotsCurr))
   intKnots[[ell]] <- intKnotsCurr
   reBlockInds[[ell]] <- currReBlockStt:ncol(Z)
   currReBlockStt <- ncol(Z) + 1
}

Ztilde <- Z*(1/sqrt(2*lambda))
yX <- as.matrix(y*X) ; yZtilde <- as.matrix(y*Ztilde)

LowRankQPobj <- LowRankQP(yZtilde,rep(-1,n),t(yX),rep(0,ncX),
                          rep(1,n),niter = 500)

# Obtain estimated beta and u vectors:

betaHat <- as.vector(LowRankQPobj$beta) 
uHat <-  as.vector(crossprod(yZtilde,LowRankQPobj$alpha))

# Plot slice of each additive model component with  
# all other variables set at their mean:

ng <- 201
xg <- seq(0,1,length = ng)
meansVec <- apply(X[,-1],2,mean)
Xmeang <- rep(1,ng)
for (j in 1:(numIndctr+numContns))
   Xmeang <- cbind(Xmeang,rep(meansVec[j],ng))

Zmeang <- NULL
for (ell in 1:numContns)
   Zmeang <- cbind(Zmeang,ZOSull(rep(meansVec[numIndctr+ell],ng),
                                     c(0,1),intKnots[[ell]]))

fitg <- vector("list",(numIndctr+numContns))
for (j in 1:(numIndctr+numContns))
{
   Xgcurr <- Xmeang
   Xgcurr[,(j + 1)] <- xg
   Zgcurr <- Zmeang
   if (j>numIndctr) 
   {
      ell <- j - numIndctr      
      Zgcurr[,reBlockInds[[ell]]] <- ZOSull(xg,c(0,1),intKnots[[ell]])
   }
   fitg[[j]] <- Xgcurr%*%betaHat + Zgcurr%*%uHat
}

# Obtain plot showing additive components of the penalized
# spline support vector machine:

par(mfrow = c(5,3),mai = c(0.25,0.25,0.1,0.03))
ylimVal <- c(-0.5,0.5)
legendCexVal <- 1.3
subInds <- sample(1:nrow(X),500,replace = FALSE)
for (j in 1:(numIndctr + numContns))
{
   xgOrigCurr <- origMins[j + 1] + xg*(origMaxs[j + 1] - origMins[j + 1])
   xlimVal <- range(xgOrigCurr)
   if (j ==15) {xlimVal <- c(2500,15000) ; xgOrigCurr[xgOrigCurr < 2500] <- NA}
   plot(0,0,type = "n",xlim = xlimVal,ylim = ylimVal,bty = "l",
        xlab = "",ylab = "")
   lines(xgOrigCurr,fitg[[j]],col = "darkgreen",lwd = 2)
   xCurr <- origMins[j + 1] + X[subInds,(j + 1)]*(origMaxs[j + 1] - origMins[j + 1])
   rug(jitter(xCurr,factor = 0.2),col = "dodgerblue",quiet = TRUE)
   abline(h = 0,col = "slateblue")
   legend("topright",legend = namesVec[j],text.col = "indianred3",cex = legendCexVal)
}

# Obtain confusion matrix and estimate the test error:

yTest <- testData[,2]
predsAllTest <-  testData[,-c(1,2)]

# Obtain X matrix in original units:

XorigTest <- cbind(rep(1,nrow(testData)),predsAllTest[,indctrInds],
                   predsAllTest[,contnsInds])
Xtest <- matrix(NA,nrow(XorigTest),ncol(XorigTest))
Xtest[,1] <- 1
for (j in 2:ncol(XorigTest))
   Xtest[,j] <- (XorigTest[,j] - origMins[j])/(origMaxs[j] - origMins[j])

Ztest <- NULL
for (ell in 1:numContns)
{
   j <- ell + numIndctr
   Ztest <- cbind(Ztest,ZOSull(Xtest[,j],c(0,1),intKnots[[ell]]))
}
ZtestTilde <- Ztest*(1/sqrt(2*lambda))

# Obtain classifications based on test data:

fHatSVM <- Xtest%*%betaHat + ZtestTilde%*%uHat
yHat <- sign(fHatSVM)

# Obtain and print out the confusion matrix and test error:

confMat <- as.matrix(table(yTest,yHat))

names(dimnames(confMat)) <- rep("",2)
dimnames(confMat)[[1]] <- c("actual good buy","actual bad buy")
dimnames(confMat)[[2]] <- c("classified good buy","classified bad buy")

cat("\n")
cat("The confusion matrix is:\n")
print(confMat)

cat("\n")
misClassError <- round(100*(sum(confMat) - sum(diag(confMat)))/sum(confMat),3)

# Obtain and estimate misclassifcation rate:

misClassRate <- (100*(sum(confMat) - sum(diag(confMat)))/sum(confMat))

cat("The estimated misclassification rate is: ",misClassRate,"%\n\n",sep = "")

############ End of carAucPenSplSVM ############

