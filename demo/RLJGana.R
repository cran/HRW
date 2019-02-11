########## R script: RLJGana ##########

# For fitting the logisitic regression with
# finite normal mixture measurement error model
# of Richardson, Leblond, Jaussent & Green (2002)
# to data from a coronary heart disease study.

# Last changed: 22 AUG 2017

# Set MCMC sample sizes:

nWarm <- 200000
nKept <- 10000
nThin <- 10

# Load in required functions and packages:

library(BRugs)
library(HRW)

# Load data:

data(CHD)
xTrue <- CHD$LDL/100
w <- CHD$TC/100
vOrig <- CHD$age
v <- vOrig/100
y <- CHD$CHD

# Divide data according to x-variable observed
# or unobserved. This is done according to the
# strategy in Richardson, Leblond, Jaussent 
# & Green (2002) in which the observed x values
# (i.e. LDL/100 values) correspond to 32 cases
# (y=1) and 40 controls (y=0) being chosen at 
# random.

set.seed(7)
n <- length(y)
caseInds <- (1:n)[y == 1]
ctrlInds <- (1:n)[y == 0]
indsObsCase <- sample(caseInds,32,replace = FALSE)
indsObsCtrl <- sample(ctrlInds,40,replace = FALSE)
indsObs <- c(indsObsCase,indsObsCtrl)
nObs <- length(indsObs)
nUnObs <- n - nObs
indsUnObs <- setdiff(1:n,indsObs)
xObs <- xTrue[indsObs]
xUnobs <- xTrue[indsUnObs]
wxObs <- w[indsObs]
wxUnObs <- w[indsUnObs]
yxObs <- y[indsObs]
yxUnObs <- y[indsUnObs]
vxObs <- v[indsObs]
vxUnObs <- v[indsUnObs]

# Fit using MCMC, starting specification of the model
# in BUGS:

RichModel <- function()
{
   for (i in 1:nObs)
   {
      logit(meanYxObs[i]) <- beta0 + beta1*xObs[i] + beta2*vxObs[i]
      meanWxObs[i] <- alpha0 + alpha1*xObs[i]
      wxObs[i] ~ dnorm(meanWxObs[i],recipSigsqW)
      yxObs[i] ~ dbern(meanYxObs[i])
   }

   for (i in 1:nUnObs)
   {
      logit(meanYxUnObs[i]) <- beta0 + beta1*xUnObs[i] + beta2*vxUnObs[i]
      meanWxUnObs[i] <- alpha0 + alpha1*xUnObs[i]
      wxUnObs[i] ~ dnorm(meanWxUnObs[i],recipSigsqW)
      yxUnObs[i] ~ dbern(meanYxUnObs[i])
   }

   for (k in 1:2)
   {
      muX[k] ~ dnorm(0.0,tauMu)
      recipSigsqX[k] ~ dgamma(0.5,recipaX[k]) 
      AxRecipSq[k] <- pow(Ax[k],-2)
      recipaX[k] ~ dgamma(0.5,AxRecipSq[k])
   }

   for (i in 1:nObs)
   {
      axObs[i] ~ dcat(omegaVec[1:2]) 
      xObs[i] ~ dnorm(muX[axObs[i]],recipSigsqX[axObs[i]])
   }

   for (i in 1:nUnObs)
   {
      axUnObs[i] ~ dcat(omegaVec[1:2]) 
      xUnObs[i] ~ dnorm(muX[axUnObs[i]],recipSigsqX[axUnObs[i]])
   }

   omega ~ dunif(0,1) ; omegaVec[1] <- omega  ; omegaVec[2] <- 1 - omega

   alpha0 ~ dnorm(0.0,tauAlpha) ; alpha1 ~ dnorm(0.0,tauAlpha)

   beta0 ~ dnorm(0.0,tauBeta) ; beta1 ~ dnorm(0.0,tauBeta)
   beta2 ~ dnorm(0.0,tauBeta)

   recipSigsqW ~ dgamma(0.5,recipaW) 
   AwRecipSq <- pow(Aw,-2) ; recipaW ~ dgamma(0.5,AwRecipSq)
   sigsqW <- 1/recipSigsqW
}
modelFileName <- file.path(getwd(),"RichModel.txt")
writeModel(RichModel,modelFileName)
   
# Do Markov chain Monte Carlo sampling via BRugs:

axObsInit <- sample(1:2,nObs,replace = TRUE)
axUnObsInit <- sample(1:2,nUnObs,replace = TRUE)

parInits <- list(list(xUnObs = rep(0,nUnObs),
                      omega = 0.5,muX = rep(0,2),recipSigsqX = rep(1,2),
                      alpha0 = 0,alpha1 = 0,beta0 = 0,beta1 = 0,beta2 = 0,
                      recipSigsqW = 1,recipaW = 1,recipaX = rep(1,2),
                      axObs = axObsInit,axUnObs = axUnObsInit))

allData <- list(nObs = nObs,nUnObs = nUnObs,vxObs = vxObs,vxUnObs = vxUnObs,
                wxObs = wxObs,wxUnObs = wxUnObs,yxObs = yxObs,yxUnObs = yxUnObs,
                xObs = xObs,Aw = 100,Ax = rep(100,2),tauAlpha = 0.01,tauBeta = 1e-10,
                tauMu = 0.01)

BRugsFit(data = allData,inits = parInits,parametersToSave
              = c("omega","muX","recipSigsqX","alpha0","alpha1",
                "beta0","beta1","beta2","xUnObs"),
                nWarm = nWarm,nKept = nKept,nThin = nThin,
                modelFile = "RichModel.txt",numChains = 1)

# Extract MCMC samples:

alpha0MCMC <- samplesSample("alpha0")
alpha1MCMC <- samplesSample("alpha1")
beta0MCMC <- samplesSample("beta0")
beta1MCMC <- samplesSample("beta1")
beta2MCMC <- samplesSample("beta2")
beta2OrigMCMC <- samplesSample("beta2")/100
omegaMCMC <- samplesSample("omega")
muX1MCMC <- samplesSample("muX[1]")
muX2MCMC <- samplesSample("muX[2]")
sigsqX1MCMC <- 1/samplesSample("recipSigsqX[1]")
sigsqX2MCMC <- 1/samplesSample("recipSigsqX[2]")

# Do summary for regression parameters:

parmsMCMC <- list(cbind(alpha0MCMC,alpha1MCMC,beta0MCMC,
                        beta1MCMC,beta2OrigMCMC))
parNamesVal <- list(expression(alpha[0]),expression(alpha[1]),expression(beta[0]),
                 expression(beta[1]),expression(beta[2]))

summMCMC(parmsMCMC,parNames = parNamesVal,columnHeadCex = 2.9)

# Plot estimated density function of partially observed predictor:

ng <- 201
mu1Hat <- mean(muX1MCMC) 
sigma1Hat <- mean(sqrt(sigsqX1MCMC))
mu2Hat <- mean(muX2MCMC) 
sigma2Hat <- mean(sqrt(sigsqX2MCMC))
xLim1 <- c(mu1Hat - 3.5*sigma1Hat,mu1Hat + 3.5*sigma1Hat)
xLow <- min(xLim1) 
xUpp <- 3

xg <- seq(xLow,xUpp,length = ng)

nMCMC <- length(muX1MCMC)

densMCMC <- matrix(0,nMCMC,ng)
for (iMCMC in 1:nMCMC)
{
   for(ig in 1:ng)
   {
      densMCMC[iMCMC,ig] <- (densMCMC[iMCMC,ig] 
            + omegaMCMC[iMCMC]*dnorm(xg[ig],muX1MCMC[iMCMC],sqrt(sigsqX1MCMC[iMCMC]))
            + (1-omegaMCMC[iMCMC])*dnorm(xg[ig],muX2MCMC[iMCMC],sqrt(sigsqX2MCMC[iMCMC]))
            )
   }
}

densBayes <- apply(densMCMC,2,mean)
densLower <- apply(densMCMC,2,quantile,0.025)
densUpper <- apply(densMCMC,2,quantile,0.975)

par(mfrow = c(1,1),mai = c(1.02,1.0,0.2,0.2))
cex.labVal <- 1.8
cex.axisVal <- 1.8
ylimVal <- range(c(densLower,densUpper))
plot(0,0,type = "n",bty = "l",xlim = range(xg),ylim = ylimVal,
     xlab = "(low density lipoprotein cholesterol level)/100",
     ylab = "estimated density function",cex.lab = cex.labVal,
     cex.axis = cex.axisVal)

polygon(c(xg,rev(xg)),c(densLower,rev(densUpper)),col = "palegreen",border = FALSE)
lines(xg,densBayes,col = "darkgreen",lwd = 2)
abline(0,0,col = "slateblue")

########### End of RLJGana ###########



