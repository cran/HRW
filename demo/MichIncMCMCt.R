########## R script: MichIncMCMC ##########

# For robust t-based nonparametric regression
# to two variables from the 1987 cross-section
# of the Michigan Panel Study of Income Dynamics.

# Last changed: 22 AUG 2017

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Load required packages:

library(rstan)
library(HRW)
library(Ecdat)

# Read in data:

data(Workinghours)
xOrig <- Workinghours$age 
yOrig <- Workinghours$income/10

# Obtain standardised data for Bayesian analysis:

mean.x <- mean(xOrig) ; mean.y <- mean(yOrig)
sd.x <- sd(xOrig) ; sd.y <- sd(yOrig)
x <- (xOrig - mean.x)/sd.x  ; y <- (yOrig - mean.y)/sd.y

# Set MCMC parameters:

nWarm <- 1000      # Length of burn-in.
nKept <- 1000        # Size of the kept sample.
nThin <- 1           # Thinning factor. 

# Set up X matrix for linear component:

X <- cbind(rep(1,length(y)),x)

# Set up Z matrix of spline basis functions:

numIntKnots <- 25
intKnots <- quantile(unique(x),seq(0,1,length = 
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])
range.x <- c(1.01*min(x) - 0.01*max(x),1.01*max(x) - 0.01*min(x))
Z <- ZOSull(x,range.x,intKnots)
ncZ <- ncol(Z)

# Specify model in Stan:

tNpRegModel <- 
   'data
   {
      int<lower = 1> n;         int<lower = 1> ncZ;
      vector[n] y;              matrix[n,2] X;
      matrix[n,ncZ] Z;          real<lower = 0> sigmaBeta;
      real<lower = 0> Au;       real<lower = 0> Aeps;
      real<lower = 0> nuLow;    real<lower = 0> nuUpp;
   }
   parameters 
   {
      vector[2] beta;            vector[ncZ] u;
      real<lower = 0> sigmaU;    real<lower = 0> sigmaEps;  
      real<lower = 0> nu;  
   }
   model 
   {
      y ~ student_t(nu,X*beta + Z*u,sigmaEps) ;
      u ~ normal(0,sigmaU);  beta ~ normal(0,sigmaBeta);
      sigmaEps ~ cauchy(0,Aeps); sigmaU ~ cauchy(0,Au);
      nu ~ uniform(nuLow,nuUpp);
   }'

allData <- list(n = length(y),ncZ = ncZ,y = y,X = X,Z = Z,
                sigmaBeta = 1e5,Au = 1e5,Aeps = 1e5,nuLow = 0.01,nuUpp = 100)

# Compile code for model if required:

if (compileCode)  
   stanCompilObj <- stan(model_code = tNpRegModel,data = allData,iter = 1,
                         chains = 1)

# Perform MCMC:

stanObj <- stan(model_code = tNpRegModel,data = allData,
                warmup = nWarm,iter = (nWarm + nKept),
                chains = 1,thin = nThin,refresh = 100,fit = stanCompilObj)

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
nuMCMC <- as.vector(extract(stanObj,"nu",permuted = FALSE))

# Display fitted curve:

cexVal <- 0.4 ; cex.labVal <- 1.8 ; lwdVal <- 2
cex.axisVal <- 1.8 
ptCol <- "dodgerblue"    ; cexVal <- 0.4
estFunCol <- "darkgreen" ; varBandCol <- "palegreen"
ng <- 101
xg <- seq(range.x[1],range.x[2],length = ng)
xgOrig <- sd.x*xg + mean.x
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,range.x,intKnots)
fMCMC <-   Xg%*%betaMCMC + Zg%*%uMCMC 
fOrigMCMC <- sd.y*fMCMC + mean.y
credLower <- apply(fOrigMCMC,1,quantile,0.025)
credUpper <- apply(fOrigMCMC,1,quantile,0.975)
fgOrig <- apply(fOrigMCMC,1,mean)

par(mfrow = c(1,2),mai = c(0.9,0.9,0.1,0.1))
ylabString <- "other household income ('000 $US)"
plot(xgOrig,fgOrig,type = "n",bty = "l",xlab = "wife's age in years",
     ylab = ylabString,xlim = range(xOrig),ylim = range(yOrig),
     cex.lab = cex.labVal,cex.axis = cex.axisVal)
points(xOrig,yOrig,col = ptCol,cex = cexVal)
polygon(c(xgOrig,rev(xgOrig)),c(credLower,rev(credUpper)),
       col = varBandCol,border = FALSE)
lines(xgOrig,fgOrig,lwd = 2,col = estFunCol)

# Do zoomed view plots:

plot(xgOrig,fgOrig,type = "n",bty = "l",xlab = "wife's age in years",
     ylab = ylabString,xlim = range(xOrig),ylim = c(0,100),
     cex.lab = cex.labVal,cex.axis = cex.axisVal)
points(xOrig,yOrig,col = ptCol,cex = cexVal)
polygon(c(xgOrig,rev(xgOrig)),c(credLower,rev(credUpper)),
       col = varBandCol,border = FALSE)
lines(xgOrig,fgOrig,lwd = 2,col = estFunCol)

# Display summary of MCMC samples for model parameters:

indQ1 <- length(xgOrig[xgOrig <= quantile(xOrig,0.25)])
fQ1MCMC <- fMCMC[indQ1,]
indQ2 <- length(xgOrig[xgOrig <= quantile(xOrig,0.50)])
fQ2MCMC <- fMCMC[indQ2,]
indQ3 <- length(xgOrig[xgOrig <= quantile(xOrig,0.75)])
fQ3MCMC <- fMCMC[indQ3,]

parms <- list(cbind(fQ1MCMC,fQ2MCMC,fQ3MCMC,sigmaEpsMCMC,nuMCMC))
parNamesVal <- list(c("est. mean funct.","at first quartile","of age"),
                    c("est. mean funct.","at second quartile","of age"),
                    c("est. mean funct.","at third quartile","of age"),
                    c(expression(sigma[epsilon])),
                    c(expression(nu)))

summMCMC(parms,parNames = parNamesVal,numerSummCex = 1.5,
         columnHeadCex = 2.9,KDEvertLine = FALSE)

############ End of MichIncMCMCt ##########
