######### R script: npRegMNAR ##########

# For performing penalized spline-based nonparametric 
# regression when the predictor is subject to missingness
# not at random.

# Last changed: 22 AUG 2017

# Set simulation setting number:
# 
# simSettNum = 1:  Missing not at random setting used in Section 4.1
#                  of Faes, Ormerod & Wand (2011),
#                  `Journal of the American Statistical Association'.
#
# simSettNum = 2: As for simSettNum = 1 but with true regression
#                 function and true error variance as in Section 
#                 6.5.2 of Harelzak, Ruppert & Wand (2017),
#                 `Semiparametric Regression with R'.

simSettNum <- 2

# Set flag for code compilation (needed if 
# running script first time in current session):

compileCode <- TRUE

# Load required packages:

library(rstan)  ;  library(HRW)

# Set MCMC parameters:

nWarm <- 10000         # Length of burn-in.
nKept <- 10000           # Size of the kept sample.
nThin <- 10              # Thinning factor. 

# Set true values of parameters and mean function:

if (simSettNum == 1)
{
   n <- 500
   fTrue <- function(x) 
      return(sin(4*pi*x^2))

   muXTrue <- 0.5 ; sigmaXTrue <- 1/6
   sigmaEpsTrue <- sqrt(0.35)
   phiTrue <- c(3,-3)
}
if (simSettNum == 2)
{
   n <- 1000
   fTrue <- function(x) 
      return(3*exp(-78*(x-0.38)^2)+exp(-200*(x-0.75)^2) - x)

   muXTrue <- 0.5 ; sigmaXTrue <- 1/6
   sigmaEpsTrue <- 0.35
   phiTrue <- c(2.95,-2.95)
}
   
# Simulate data with missingess not at random
# according to a linear logistic regression model:

set.seed(2)
x <- rnorm(n,muXTrue,sigmaXTrue)

# Obtain missingness indicators:

probxObs <- 1/(1+exp(-(phiTrue[1]+phiTrue[2]*x)))
indicxObsTmp <- rbinom(n,1,probxObs)

# Obtain observed data:

xObs <- x[indicxObsTmp == 1]
nObs <- length(xObs)
yxObs <- fTrue(xObs) + sigmaEpsTrue*rnorm(nObs)
indicxObs <- rep(1,nObs)

# Obtain missing data:

xUnobsTrue <- x[indicxObsTmp == 0]
nUnobs <- length(xUnobsTrue)
yxUnobs <- fTrue(xUnobsTrue) +  sigmaEpsTrue*rnorm(nUnobs)
indicxUnobs <- rep(0,nUnobs)

r <- c(indicxObs,indicxUnobs)

# Obtain the Z matrix for the observed data:

ncZ <- 30
knots <- seq(min(xObs),max(xObs),length = (ncZ+2))[-c(1,ncZ+2)]
ZxObs <- outer(xObs,knots,"-")
ZxObs <- ZxObs*(ZxObs>0)

# Specify model in Stan:

npRegMNARModel <- 
   'data
   {
      int<lower=1> nObs;        int<lower=1> nUnobs;
      int<lower=1> n;           int<lower=1> ncZ;
      vector[nObs] yxObs;       vector[nUnobs] yxUnobs;
      vector[nObs] xObs;
      vector[ncZ]  knots;       matrix[nObs,ncZ] ZxObs;
      int<lower=0,upper=1> r[n]; 
      real<lower=0> sigmaBeta;  real<lower=0> sigmaMu;
      real<lower=0> sigmaPhi;
      real<lower=0> Ax;         real<lower=0> Aeps;
      real<lower=0> Au;         
   }
   transformed data
   {
      vector[n] y;
      for (i in 1:nObs)
         y[i] = yxObs[i]; 
      for (i in 1:nUnobs)
         y[i+nObs] = yxUnobs[i];
   }
   parameters 
   {
      vector[2] beta;           vector[ncZ] u;
      vector[2] phi;
      real muX;                 real<lower=0> sigmaX;
      real<lower=0> sigmaEps;   real<lower=0> sigmaU;
      real xUnobs[nUnobs];
   }
   transformed parameters 
   {
      matrix[n,2] X;       matrix[n,ncZ] Z;
      for (i in 1:nObs)
      {
         X[i,1] = 1   ;   X[i,2] = xObs[i] ;  
         Z[i] = ZxObs[i];
      }
      for (i in 1:nUnobs)
      {
         X[i+nObs,1] = 1    ;   X[i+nObs,2] = xUnobs[i];
         for (k in 1:ncZ)   
            Z[i+nObs,k] = (xUnobs[i]-knots[k])*step(xUnobs[i]-knots[k]);
      }      
   }
   model 
   {
      y ~ normal(X*beta+Z*u,sigmaEps); 
      r ~ bernoulli_logit(X*phi);
      col(X,2) ~ normal(muX,sigmaX); 
      u ~ normal(0,sigmaU) ; beta ~ normal(0,sigmaBeta);
      muX ~ normal(0,sigmaMu); phi ~ normal(0,sigmaPhi);
      sigmaX ~ cauchy(0,Ax);
      sigmaEps ~ cauchy(0,Aeps);   sigmaU ~ cauchy(0,Au);
   }'

# Set up input data:

allData <- list(nObs = nObs,nUnobs = nUnobs,n = (nObs+nUnobs),
                ncZ = ncZ,xObs = xObs,yxObs = yxObs,yxUnobs = yxUnobs,knots = knots,
                ZxObs = ZxObs,r = r,sigmaMu = 1e5,sigmaBeta = 1e5,sigmaPhi = 1e5,
                sigmaEps = 1e5,sigmaX = 1e5,Ax = 1e5,Aeps = 1e5,Au = 1e5)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = npRegMNARModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = npRegMNARModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 25,
                 fit = stanCompilObj)

# Extract relevant MCMC sampes:

betaMCMC <- NULL
for (j in 1:2)
{
   charVar <- paste("beta[",as.character(j),"]",sep  = "") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar,permuted = FALSE))
}
uMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("u[",as.character(k),"]",sep = "") 
   uMCMC <- rbind(uMCMC,extract(stanObj,charVar,permuted = FALSE))
}
muXMCMC <- as.vector(extract(stanObj,"muX",permuted = FALSE))
sigmaXMCMC <- as.vector(extract(stanObj,"sigmaX",permuted = FALSE))
sigmaEpsMCMC <- as.vector(extract(stanObj,"sigmaEps",permuted = FALSE))
sigmaUMCMC <- as.vector(extract(stanObj,"sigmaU",permuted = FALSE))
phiMCMC <- NULL
for (j in 1:2)
{
   charVar <- paste("phi[",as.character(j),"]",sep = "") 
   phiMCMC <- rbind(phiMCMC,extract(stanObj,charVar,permuted = FALSE))
}
xUnobsMCMC <- NULL
for (i in 1:nUnobs)
{
   charVar <- paste("xUnobs[",as.character(i),"]",sep = "") 
   xUnobsMCMC <- rbind(xUnobsMCMC,extract(stanObj,charVar,permuted = FALSE))
}

# Plot fit:

cexVal <- 0.3
obsCol <- "darkblue" ; misCol <- "lightskyblue"
estFunCol <- "darkgreen"; trueFunCol <- "indianred3"
varBandCol <- "palegreen"
ng <- 101
xg <- seq(min(x),max(x),length = ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- outer(xg,knots,"-")
Zg <- Zg*(Zg>0)
fMCMC <-   Xg%*%betaMCMC + Zg%*%uMCMC 
credLower <- apply(fMCMC,1,quantile,0.025)
credUpper <- apply(fMCMC,1,quantile,0.975)
fg <- apply(fMCMC,1,mean)

par(mfrow = c(1,1),mai = c(1.02,0.82,0.82,0.42))
plot(xg,fg,type = "n",bty = "l",xlab = "x",ylab = "y",xlim = range(c(xObs,xUnobsTrue)),
     ylim = range(c(yxObs,yxUnobs)))
polygon(c(xg,rev(xg)),c(credLower,rev(credUpper)),col = varBandCol,border = FALSE)

lines(xg,fg,lwd = 2,col = estFunCol)
lines(xg,fTrue(xg),lwd = 2,col = trueFunCol)
points(xObs,yxObs,col = obsCol,cex = cexVal)
points(xUnobsTrue,yxUnobs,col = misCol,cex = cexVal)
legend(0.75,1.5,c("fully observed","x value unobserved"),col = c(obsCol,misCol),
       pch = rep(1,2),pt.lwd = rep(2,2))
legend(0.467,-1.2,c("true f","estimated f"),col = c(trueFunCol,estFunCol),
       lwd = rep(2,2))

# Display summary of MCMC samples for model parameters:

indQ1 <- length(xg[xg <= quantile(xObs,0.25)])
fQ1MCMC <- fMCMC[indQ1,]
indQ2 <- length(xg[xg <= quantile(xObs,0.50)])
fQ2MCMC <- fMCMC[indQ2,]
indQ3 <- length(xg[xg <= quantile(xObs,0.75)])
fQ3MCMC <- fMCMC[indQ3,]

parms <- list(cbind(muXMCMC,sigmaXMCMC,sigmaEpsMCMC,phiMCMC[1,],phiMCMC[2,],
                    fQ1MCMC,fQ2MCMC,fQ3MCMC))
parNamesVal <- list(c(expression(mu[x])),c(expression(sigma[x])),
                      c(expression(sigma[epsilon])),
                      c(expression(phi[0])),c(expression(phi[1])),
                      c("first quartile","of x"),
                      c("second quart.","of x"),
                      c("third quartile","of x"))

summMCMC(parms,parNames = parNamesVal,numerSummCex = 1,
            KDEvertLine = FALSE,addTruthToKDE = c(muXTrue,sigmaXTrue,sigmaEpsTrue,
            phiTrue[1],phiTrue[2],fTrue(xg[indQ1]),
            fTrue(xg[indQ2]),fTrue(xg[indQ3])))

# Display summary for the first five missing x values:

parms <- list(t(xUnobsMCMC[1:5,]))
parNamesVal <- list(c(expression(x[1]^{"unobs"})),
                      c(expression(x[2]^{"unobs"})),
                      c(expression(x[3]^{"unobs"})),
                      c(expression(x[4]^{"unobs"})),
                      c(expression(x[5]^{"unobs"})))

summMCMC(parms,parNames = parNamesVal,KDEvertLine = FALSE,addTruthToKDE = x[1:5])

########## End of npRegMNAR ##########
