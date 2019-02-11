######### R script: npRegMCAR ##########

# For performing penalised spline-based nonparametric 
# regression when the predictor is subject to missingness
# completely at random.

# Last changed: 22 AUG 2017

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Load required package:

library(rstan) ; library(HRW)

# Set MCMC parameters:

nWarm <- 1000        # Length of burn-in.
nKept <- 1000          # Size of the kept sample.
nThin <- 1             # Thinning factor. 

# Set true values of mean function and model parameters:

fTrue <- function(x) 
   return(3*exp(-78*(x-0.38)^2)+exp(-200*(x-0.75)^2)-x)

muXTrue <- 0.5 ; sigmaXTrue <- 1/6
sigmaEpsTrue <- 0.35

# Simulate data using these values; with
# 20% of x data missing completely at random:

set.seed(1)
n <- 1000
misPropn <- 0.2
x <- rnorm(n,muXTrue,sigmaXTrue)
y <- fTrue(x) + sigmaEpsTrue*rnorm(n)
nUnobs <- round(misPropn*n) ; nObs <- n - nUnobs
indsUnobs <- sample(1:n,nUnobs,replace=FALSE)
indsObs <- setdiff(1:n,indsUnobs)
xUnobsTrue <- x[indsUnobs] ; xObs <- x[indsObs]
yxUnobs <- y[indsUnobs] ; yxObs <- y[indsObs]

# Obtain the Z matrix for the observed data:

ncZ <- 25
knots <- seq(min(xObs),max(xObs),length=(ncZ+2))[-c(1,ncZ+2)]
ZxObs <- outer(xObs,knots,"-")
ZxObs <- ZxObs*(ZxObs>0)

# Specify model in Stan:

npRegMCARModel <- 
   'data
   {
      int<lower=1> nObs;       int<lower=1> nUnobs;
      int<lower=1> ncZ;
      vector[nObs] yxObs;       vector[nUnobs] yxUnobs;
      vector[nObs] xObs;
      vector[ncZ]  knots;       matrix[nObs,ncZ] ZxObs;
      real<lower=0> sigmaBeta;  real<lower=0> sigmaMu;
      real<lower=0> Ax;         real<lower=0> Aeps;
      real<lower=0> Au;         
   }
   parameters 
   {
      vector[2] beta;           vector[ncZ] u;
      real muX;                 real<lower=0> sigmaX;
      real<lower=0> sigmaEps;   real<lower=0> sigmaU;
      real xUnobs[nUnobs];
   }
   transformed parameters 
   {
      matrix[nObs,2] XxObs;    matrix[nUnobs,2] XxUnobs;
      matrix[nUnobs,ncZ] ZxUnobs;
      for (i in 1:nObs)
      {
         XxObs[i,1] = 1    ;   XxObs[i,2] = xObs[i];
      }
      for (i in 1:nUnobs)
      {
         XxUnobs[i,1] = 1    ;   XxUnobs[i,2] = xUnobs[i];
         for (k in 1:ncZ)   
         {
            ZxUnobs[i,k] = (xUnobs[i]-knots[k])*step(xUnobs[i]-knots[k]);
         }
      }
   }
   model 
   {
      yxObs ~ normal(XxObs*beta+ZxObs*u,sigmaEps);
      xObs  ~ normal(muX,sigmaX); 
      yxUnobs ~ normal(XxUnobs*beta+ZxUnobs*u,sigmaEps); 
      xUnobs  ~ normal(muX,sigmaX); 
      u ~ normal(0,sigmaU) ; beta ~ normal(0,sigmaBeta);
      muX ~ normal(0,sigmaMu); sigmaX ~ cauchy(0,Ax);
      sigmaEps ~ cauchy(0,Aeps);   sigmaU ~ cauchy(0,Au);
   }'

# Set up input data:

allData <- list(nObs = nObs,nUnobs = nUnobs,ncZ = ncZ,xObs = xObs,
                yxObs = yxObs,yxUnobs = yxUnobs,ZxObs = ZxObs,knots = knots,
                sigmaMu = 1e5,sigmaBeta = 1e5,sigmaEps = 1e5,sigmaX = 1e5,Ax = 1e5,
                Aeps = 1e5,Au = 1e5)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = npRegMCARModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = npRegMCARModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 100,
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
xUnobsMCMC <- NULL
for (i in 1:nUnobs)
{
   charVar <- paste("xUnobs[",as.character(i),"]",sep = "") 
   xUnobsMCMC <- rbind(xUnobsMCMC,extract(stanObj,charVar,permuted = FALSE))
}

# Plot fit:

obsCol <- "darkblue" ; misCol <- "lightskyblue"
estFunCol <- "darkgreen"; trueFunCol <- "indianred3"
varBandCol <- "palegreen" ; cexVal <- 0.3
ng <- 201
xg <- seq(0,1,,length = ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- outer(xg,knots,"-")
Zg <- Zg*(Zg>0)
fMCMC <-   Xg%*%betaMCMC + Zg%*%uMCMC 
credLower <- apply(fMCMC,1,quantile,0.025)
credUpper <- apply(fMCMC,1,quantile,0.975)
fg <- apply(fMCMC,1,mean)

par(mfrow = c(1,1),mai = c(1.02,0.82,0.82,0.42))
plot(xg,fg,type = "n",bty = "l",xlab = "x",ylab = "y",xlim = range(x),ylim = range(y))
polygon(c(xg,rev(xg)),c(credLower,rev(credUpper)),col = varBandCol,border = FALSE)

lines(xg,fg,lwd = 2,col = estFunCol)
lines(xg,fTrue(xg),lwd = 2,col = trueFunCol)
points(xObs,yxObs,col = obsCol,cex = cexVal)
points(xUnobsTrue,yxUnobs,col = misCol,cex = cexVal)
legend("topright",c("fully observed","x value unobserved"),col = c(obsCol,misCol),
       pch = rep(1,2),pt.cex = cexVal)
legend("topleft",c("true f","estimated f"),col = c(trueFunCol,estFunCol),
       lwd = rep(2,2))

# Display summary of MCMC samples for model parameters:

indQ1 <- length(xg[xg <= quantile(xObs,0.25)])
fQ1MCMC <- fMCMC[indQ1,]
indQ2 <- length(xg[xg <= quantile(xObs,0.50)])
fQ2MCMC <- fMCMC[indQ2,]
indQ3 <- length(xg[xg <= quantile(xObs,0.75)])
fQ3MCMC <- fMCMC[indQ3,]

parms <- list(cbind(muXMCMC,sigmaXMCMC,sigmaEpsMCMC,fQ1MCMC,fQ2MCMC,fQ3MCMC))
parNamesVal <- list(c(expression(mu[x])),c(expression(sigma[x])),
                      c(expression(sigma[epsilon])),
                      c("first quartile","of x"),
                      c("second quart.","of x"),
                      c("third quartile","of x"))

summMCMC(parms,parNames = parNamesVal,numerSummCex = 1,
         KDEvertLine = FALSE,addTruthToKDE = c(muXTrue,sigmaXTrue,
         sigmaEpsTrue,fTrue(xg[indQ1]),
         fTrue(xg[indQ2]),fTrue(xg[indQ3])))

# Display summary for some of the missing x values:

parms <- list(t(xUnobsMCMC[1:5,]))
parNamesVal <- list(c(expression(x["unobs,1"])),
                      c(expression(x["unobs,2"])),
                      c(expression(x["unobs,3"])),
                      c(expression(x["unobs,4"])),
                      c(expression(x["unobs,5"])))
summMCMC(parms,parNames = parNamesVal,KDEvertLine = FALSE,addTruthToKDE = x[1:5])

############ End of npRegMCAR ############
