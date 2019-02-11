######### R script: MichIncMCMCskewt ##########

# For fitting a Jones-Faddy skew-t distribution
# via MCMC to two variables from the 1987 cross-section
# of the Michigan Panel Study of Income Dynamics.

# Last changed: 22 AUG 2017

# Set flag for code compilation (needed if 
# running script first time in current session):

compileCode <- TRUE

# Load required packages:

library(HRW) ; library(rstan) ; 
library(Ecdat)

# Read in data:

data(Workinghours)
xOrig <- Workinghours$age 
yOrig <- Workinghours$income/10
n <- length(yOrig)

# Obtain standardised data for Bayesian analysis:

mean.x <- mean(xOrig) ; mean.y <- mean(yOrig)
sd.x <- sd(xOrig) ; sd.y <- sd(yOrig)
x <- (xOrig - mean.x)/sd.x  ; y <- (yOrig - mean.y)/sd.y

# Set MCMC parameters:

nWarm <- 1000        # Length of burn-in.
nKept <- 1000          # Size of the kept sample.
nThin <- 1             # Thinning factor. 

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

JonesFaddyModel <- 
   'data
   {
      int<lower = 1> n;       int<lower = 1> ncZ;
      vector[n] y;            matrix[n,2] X;
      matrix[n,ncZ] Z;        real<lower = 0> sigmaBeta;
      real<lower = 0> Au;     real<lower = 0> Aeps;
   }
   parameters 
   {
      vector[2] beta;            vector[ncZ] u;             
      real<lower = 0> sigmaEps;  real<lower = 0> sigmaU;
      real<lower = 0> nuLft;     real<lower = 0> nuRgt;
   }
   transformed parameters
   { 
      real JFmode;           real logNormFac;
      vector[n] arg;
      JFmode = (nuLft-nuRgt)*sqrt((nuLft + nuRgt)/((2*nuLft + 1)*(2*nuRgt + 1)));
      arg = JFmode + (y - X*beta - Z*u)/sigmaEps;   
      logNormFac = (nuLft+nuRgt-1)*log(2) + 0.5*log(nuLft + nuRgt)
                    - lgamma(nuLft+nuRgt) + lgamma(nuLft) + lgamma(nuRgt);
   }
   model 
   {
      for (i in 1:n)
         target += (nuLft+0.5)*log1p(arg[i]/sqrt(nuLft+nuRgt+square(arg[i])))
                         + (nuRgt+0.5)*log1m(arg[i]/sqrt(nuLft+nuRgt+square(arg[i])))
                         - log(sigmaEps) - logNormFac;
      u ~ normal(0,sigmaU); beta ~ normal(0,sigmaBeta);   
      sigmaU ~ cauchy(0,Au); sigmaEps ~ cauchy(0,Aeps);
      nuLft ~ uniform(0.01,100); nuRgt ~ uniform(0.01,100);
   }'

# Store data in a list in format required by Stan:

X <- cbind(rep(1,length(y)),x)
allData <- list(n = n,ncZ = ncZ,y = y,X = X,Z = Z,sigmaBeta = 1e5,Aeps = 1e5,Au = 1e5)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = JonesFaddyModel,data = allData,
                         iter = 1,chains = 1)

# Obtain MCMC samples for each parameter using Stan:

stanObj <-  stan(model_code = meanVarModel,data = allData,
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
nuLftMCMC <- as.vector(extract(stanObj,"nuLft",permuted = FALSE))
nuRgtMCMC <- as.vector(extract(stanObj,"nuRgt",permuted = FALSE))

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

parms <- list(cbind(fQ1MCMC,fQ2MCMC,fQ3MCMC,sigmaEpsMCMC,nuLftMCMC,nuRgtMCMC))
parNamesVal <- list(c("mode function","at first quartile","of age"),
                    c("mode function","at second quartile","of age"),
                    c("mode function","at third quartile","of age"),
                    c(expression(sigma[epsilon])),
                    c(expression(nu[left])),c(expression(nu[right])))

summMCMC(parms,parNames = parNamesVal,numerSummCex = 1.5,KDEvertLine = FALSE)

########## End of MichIncMCMCskewt ##########
