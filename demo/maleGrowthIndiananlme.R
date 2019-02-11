########## R script: maleGrowthIndiananlme ##########

# For performing a frequentist group-specific curve analysis
# of the Indiana growth data for male subjects using the
# R function lme() within the package nlme.

# Last changed: 19 JUN 2017

# Load required packages:

library(HRW)  ; library(lattice)  ;  library(nlme)  

# Load in the data:

data(growthIndiana)

# Extract data on males:

growthIndianaMales <- growthIndiana[growthIndiana$male==1,]

# Create relevant variables: 

yOrig <- growthIndianaMales$height
xOrig <- growthIndianaMales$age
idnumOrig <- growthIndianaMales$idnum
typeIsB <- growthIndianaMales$black

# Create new (ordered) ID numbers:

idnum <- rep(NA,length(idnumOrig))
uqID <- unique(idnumOrig)
for (i in 1:length(uqID))
      idnum[idnumOrig==uqID[i]] <- i

# Save mean and standard deviation information:

mean.x <- mean(xOrig)   ;   sd.x <- sd(xOrig)
mean.y <- mean(yOrig)   ;   sd.y <- sd(yOrig)

# Store the standardised versions of the predictor
# and response variables in x and y:

x <- (xOrig - mean.x)/sd.x
y <- (yOrig - mean.y)/sd.y

# Do plot of raw data:

allDataDF <- data.frame(xOrig,yOrig,idnum,typeIsB)
cex.labVal <- 1.35
colourVec <- c("dodgerblue","deeppink","blue","darkmagenta")

figRaw <- xyplot(yOrig~xOrig|idnum,groups=idnum,data=allDataDF,
                 layout=c(9,13),xlab=list(label="age (years)",cex=cex.labVal),
                 ylab=list(label="height (centimetres)",cex=cex.labVal),
                 scales=list(cex=1.25),strip=FALSE,as.table=TRUE,
                 key = list(title="",
                            columns = 2, 
                            points = list(pch = rep(1,2),col = colourVec[1:2]),
                            text = list(c("black","white"),cex=1.55)),
                 panel=function(x,y,subscripts,groups)
                 {
                    panel.grid()
                    colourInd <- 2 - typeIsB[subscripts[1]] 
                    panel.superpose(x,y,subscripts,groups,
                             type="b",col=colourVec[colourInd],pch=1,cex=0.5)

                 })
print(figRaw)

# Set dimension variables:

numObs <- length(x)
numGrp <- length(unique(idnum))

# Set up X and Z matrix for mean curve:

Xbase <- cbind(rep(1,numObs),x)
XBvA <- typeIsB*Xbase
X <- cbind(Xbase,XBvA) 

numIntKnots <- 20
intKnots <- quantile(unique(x),
                     seq(0,1,length=numIntKnots+2))[-c(1,numIntKnots+2)]

range.x <- c(1.01*min(x) - 0.01*max(x),1.01*max(x)-0.01*min(x))

Zbase <- ZOSull(x,intKnots=intKnots,range.x=range.x)
ZBvA <- typeIsB*Zbase

# Set up Z matrix for group specific curves:

numIntKnotsGrp <- 10
intKnotsGrp <- quantile(unique(x),
        seq(0,1,length=numIntKnotsGrp+2))[-c(1,numIntKnotsGrp+2)]

Zgrp <- ZOSull(x,intKnots=intKnotsGrp,range.x=range.x)

# Let up dimension variables:

ncZ <- ncol(Zbase) ; ncZgrp <- ncol(Zgrp)

# Fit using linear mixed model software:

dummyId <- factor(rep(1,numObs))
Zblock <- list(dummyId=pdIdent(~-1+Zbase),
               dummyId=pdIdent(~-1+ZBvA),
               idnum=pdSymm(~x),
               idnum=pdIdent(~-1+Zgrp))
growthINmalGD <- groupedData(y~X[,-1]|rep(1,length=numObs),
                         data=data.frame(y,X,idnum))
cat("Starting fitting (takes several minutes)...\n")
fit <- lme(y~-1+X,data=growthINmalGD,random=Zblock)
cat("Finished fitting.\n")

sig.epsHat <- fit$sigma
sig.uHat <- sig.epsHat*exp(unlist(fit$modelStruct))

lamVal <- (sig.epsHat/sig.uHat)^2

# Set up plotting grids:

ng <- 101
xg <- seq(range.x[1],range.x[2],length=ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,intKnots=intKnots,range.x=range.x)
betaHat <- as.vector(fit$coef$fixed)
uHatBase <- as.vector(fit$coef$random[[1]])
uHatCont <- as.vector(fit$coef$random[[2]])

fhatBaseg <- Xg%*%betaHat[1:2] + Zg%*%uHatBase
Contg <-  Xg%*%betaHat[3:4] + Zg%*%uHatCont
fhatBlackg <- fhatBaseg + Contg

# Plot fitted curves:

Xsubjg <- cbind(rep(1,ng),xg)
Zgrpg <- ZOSull(xg,intKnots=intKnotsGrp,range.x=range.x)
meanCrvs <- vector("list",numGrp)
for (i in 1:numGrp)
{
   indsi <- (1:numObs)[idnum==i]
   uLinHati <- as.vector(fit$coef$random[[3]][i,])
   uSplHati <- as.vector(fit$coef$random[[4]][i,])
   ghati <- Xsubjg%*%uLinHati + Zgrpg%*%uSplHati
   meanCrvs[[i]] <- fhatBaseg + typeIsB[indsi[1]]*Contg + ghati
}

# Convert to original units:

xgOrig <- xgOrig <- mean.x + sd.x*xg
meanCrvsOrig <- vector("list",numGrp)
for (i in 1:numGrp)
   meanCrvsOrig[[i]] <- mean.y + sd.y*meanCrvs[[i]]
   
# Do lattice type plot showing group-specific fits:

figFit <- xyplot(yOrig~xOrig|idnum,groups=idnum,data=allDataDF,
                 layout=c(9,13),xlab=list(label="age (years)",cex=cex.labVal),
                 ylab=list(label="height (centimetres)",cex=cex.labVal),
                 scales=list(cex=1.25),strip=FALSE,as.table=TRUE,
                 key = list(title="",
                            columns = 2, 
                            lines = list(lty = rep(1,2),col = colourVec[3:4]),
                            text = list(c("black","white"),cex=1.55)),
                 panel=function(x,y,subscripts,groups)
                 {
                    panel.grid() 
                    iGrp <- idnum[subscripts][1]
                    colourInd <- 2 - typeIsB[subscripts[1]]
                    panel.superpose(x,y,subscripts,groups,
                             type="b",col=colourVec[colourInd],pch=1,cex=0.5)
                    panel.xyplot(xgOrig,meanCrvsOrig[[iGrp]],
                                 type="l",col=colourVec[colourInd+2])
                 })
print(figFit)

# Obtain penalty matrix:

sigsq.epsHat <- fit$sigma^2
sig.uHat <- as.numeric(sqrt(sigsq.epsHat)*exp(unlist(fit$modelStruct)))

sigsq.Gbl <- sig.uHat[6]^2
sigsq.Gblcont <- sig.uHat[5]^2
cholSigma <- rbind(c(sig.uHat[2],sig.uHat[3]),c(0,sig.uHat[4]))
SigmaHat <- crossprod(cholSigma)
sigsq.Sbj <- sig.uHat[1]^2

DmatGbl <- (sigsq.epsHat/sigsq.Gbl)*diag(ncZ)
DmatGblcont <- (sigsq.epsHat/sigsq.Gblcont)*diag(ncZ)
DmatLinSbj <- sigsq.epsHat*kronecker(diag(numGrp),solve(SigmaHat))
DmatSplSbj <- (sigsq.epsHat/sigsq.Sbj)*diag(ncol(Zgrp)*numGrp)

dimVec <- c(4,nrow(DmatGbl),nrow(DmatGblcont),nrow(DmatLinSbj),
             nrow(DmatSplSbj))

lamMat <- matrix(0,sum(dimVec),sum(dimVec))
csdV <- cumsum(dimVec)
lamMat[(csdV[1]+1):csdV[2],(csdV[1]+1):csdV[2]] <- DmatGbl
lamMat[(csdV[2]+1):csdV[3],
       (csdV[2]+1):csdV[3]] <- DmatGblcont
lamMat[(csdV[3]+1):csdV[4],
       (csdV[3]+1):csdV[4]] <- DmatLinSbj
lamMat[(csdV[4]+1):csdV[5],
       (csdV[4]+1):csdV[5]] <- DmatSplSbj

# Obtain C matrix:

uqID <- unique(idnum)
Cmat <- cbind(X,Zbase,ZBvA)
for (iSbj in 1:numGrp)
{
   newCols <- matrix(0,numObs,2)
   indsCurr <- (1:numObs)[idnum==uqID[iSbj]]
   newCols[indsCurr,] <- X[indsCurr,1:2]
   Cmat <- cbind(Cmat,newCols)
}
for (iSbj in 1:numGrp)
{
   newCols <- matrix(0,numObs,ncol(Zgrp))
   indsCurr <- (1:numObs)[idnum==uqID[iSbj]]
   newCols[indsCurr,] <- Zgrp[indsCurr,]
   Cmat <- cbind(Cmat,newCols)
}

# Obtain full covariance matrix:

CTC <- crossprod(Cmat)
fullCovMat <- solve(CTC + lamMat)

# Find subset of covariance corresponding to the contrast curve:

contInds <- c(3,4,(ncZ+5):(2*ncZ+4))
contCovMat <- fullCovMat[contInds,contInds]

# Obtain approximate pointwise 95% confidence limits:

Cg <- cbind(Xg,Zg)
sdg <- sqrt(sigsq.epsHat)*sqrt(diag(Cg%*%contCovMat%*%t(Cg)))
lowerg <- Contg - qnorm(0.975)*sdg ; upperg <- Contg + qnorm(0.975)*sdg

ContgOrig <- sd.y*Contg
lowergOrig <- sd.y*lowerg
uppergOrig <- sd.y*upperg

# Plot contrast with variability band:

par(mai=c(1,0.9,0.3,0.2))
cex.labVal <- 2   ; cex.axisVal <- 1.5
plot(xgOrig,ContgOrig,bty="l",
     type="n",xlab="age (years)",ylab="mean difference in height (cm)",
     ylim=c(-5,10),cex.lab=cex.labVal,cex.axis=cex.axisVal)

polygon(c(xgOrig,rev(xgOrig)),c(lowergOrig,rev(uppergOrig)),col="PaleGreen",
     border=FALSE)

lines(xgOrig,ContgOrig,col="DarkGreen",lwd=2)
abline(0,0,col="slateblue")
rug(xOrig,col="dodgerblue")

############ End of maleGrowthIndiananlme ###########








