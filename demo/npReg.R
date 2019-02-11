########## R script: npReg ##########

# For performing MCMC-based penalized spline
# nonparametric regression using Stan.

# Last changed: 22 AUG 2017

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Load required packages:

library(HRW) ; library(rstan) ; library(KernSmooth)

# Set MCMC parameters:

nWarm <- 1000      # Length of burn-in.
nKept <- 1000        # Size of the kept sample.
nThin <- 1           # Thinning factor. 

# Set true values of mean function and model parameters:

fTrue <- function(x) 
   return(3*exp(-78*(x-0.38)^2)+exp(-200*(x-0.75)^2) - x)
muXTrue <- 0.5 ; sigmaXTrue <- 1/6
sigmaEpsTrue <- 0.35

# Simulate data according to classic measurement error
# model:

set.seed(2)
n <- 1000
x <- rnorm(n,muXTrue,sigmaXTrue)
y <- fTrue(x) + sigmaEpsTrue*rnorm(n)

# Obtain knots and spline basis functions for x variable:

a <- 1.01*min(x) - 0.01*max(x)  
b <- 1.01*max(x) - 0.01*min(x)  
numIntKnots <- 25
intKnots <-  quantile(unique(x),seq(0,1,length = numIntKnots+2)
                      [-c(1,numIntKnots+2)])
Z <- ZOSull(x,intKnots = intKnots,range.x = c(a,b))
ncZ <- ncol(Z)

# Specify model in Stan:

npRegModel <- 
   'data
   {
      int<lower = 1> n;         int<lower = 1> ncZ;
      vector[n] y;            matrix[n,2] X;
      matrix[n,ncZ] Z;        real<lower = 0> sigmaBeta;
      real<lower = 0> Au;       real<lower = 0> Aeps;
   }
   parameters 
   {
      vector[2] beta;          vector[ncZ] u;
      real<lower = 0> sigmaEps;  real<lower = 0> sigmaU;
   }
   model 
   {
      y ~ normal(X*beta + Z*u,sigmaEps);
      u ~ normal(0,sigmaU); beta ~ normal(0,sigmaBeta); 
      sigmaEps ~ cauchy(0,Aeps); sigmaU ~ cauchy(0,Au);
   }'

# Store data in a list in format required by Stan:

X <- cbind(rep(1,length(y)),x)
allData <- list(n = length(x),ncZ = ncZ,y = y,X = X,Z = Z,
                sigmaBeta = 1e5,Au = 1e5,Aeps = 1e5)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = npRegModel,data = allData,
                         iter = 1,chains = 1)

# Obtain MCMC samples for each parameter using Stan:

stanObj <-  stan(model_code = npRegModel,data = allData,
                 warmup = nWarm,iter = (nWarm + nKept),
                 chains = 1,thin = nThin,refresh = 100,
                 fit = stanCompilObj)

# Extract relevant MCMC samples:

betaMCMC <- NULL
for (j in 1:2)
{
   charVar <- paste("beta[",as.character(j),"]",sep = "") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar,permuted = FALSE))
}
uMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("u[",as.character(k),"]",sep = "") 
   uMCMC <- rbind(uMCMC,extract(stanObj,charVar,permuted = FALSE))
}
sigmaEpsMCMC <- as.vector(extract(stanObj,"sigmaEps",permuted = FALSE))
sigmaUMCMC <- as.vector(extract(stanObj,"sigmaU",permuted = FALSE))

# Obtain MCMC samples of regression curves over a fine grid:

ng <- 201
xg <- seq(a,b,length = ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,intKnots = intKnots,range.x = c(a,b))
fhatMCMC <- Xg%*%betaMCMC + Zg%*%uMCMC

# Convert fhatMCMC matrix to original scale:

fhatg <- apply(fhatMCMC,1,mean)
credLower <- apply(fhatMCMC,1,quantile,0.025)
credUpper <- apply(fhatMCMC,1,quantile,0.975)

# Display the fit:

par(mfrow = c(1,1))

cexVal <- 0.3
estFunCol <- "darkgreen"; trueFunCol <- "indianred3"
varBandCol <- "palegreen"
plot(x,y,type = "n",xlab = "x",ylab = "y",bty = "l")
polygon(c(xg,rev(xg)),c(credLower,rev(credUpper)),
        col = varBandCol,border = FALSE)
lines(xg,fhatg,col = estFunCol,lwd = 2)
points(x,y,col = "darkblue",cex = cexVal)
lines(xg,fTrue(xg),lwd = 2,col = trueFunCol)

# Do some summaries and diagnostic checking of the MCMC:

indQ1 <- length(xg[xg<quantile(x,0.25)])
indQ2 <- length(xg[xg<quantile(x,0.50)])
indQ3 <- length(xg[xg<quantile(x,0.75)])
fhatQ1 <- fhatMCMC[indQ1,]
fhatQ2 <- fhatMCMC[indQ2,]
fhatQ3 <- fhatMCMC[indQ3,]
parms <- list(cbind(sigmaEpsMCMC,fhatQ1,fhatQ2,fhatQ3))
parNamesVal <- list(c(expression(sigma[epsilon])),
                      c("funct. est. at 1st","quantile of",
                        "predictor"),
                      c("funct. est. at 2nd","quantile of",
                        "predictor"),
                      c("funct. est. at 3rd","quantile of",
                        "predictor"))

summMCMC(parms,parNames = parNamesVal,numerSummCex = 1.25,
         columnHeadCex = 2.9,KDEvertLine = FALSE,
         addTruthToKDE = c(sigmaEpsTrue,fTrue(xg[indQ1]),
         fTrue(xg[indQ2]),fTrue(xg[indQ3])))

########## End of npReg ##########

