######### R script: npRegBerksonMeaErr ##########

# For performing penalized spline-based nonparametric 
# regression when the predictor is subject to classical
# measurement error.

# Last changed: 28 NOV 2017

cat("\n
    WARNING: tests for this example indicate very high\n
    Markov chain Monte Carlo sample sizes are required\n
    for good chain diagnostic plots. Therefore, this\n
    script may take very many hours to run.\n\n\n")

readline("Hit Enter to continue.\n")

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Load required packages:

library(rstan)  ;  library(HRW)

# Set MCMC parameters:

nWarm <- 50000             # Length of burn-in.
nKept <- 10000             # Size of the kept sample.
nThin <- 10                # Thinning factor. 

# Set true values of mean function and model parameters:

fTrue <- function(x) 
   return(3*exp(-78*(x-0.38)^2)+exp(-200*(x-0.75)^2)-x)
  
muWTrue <- 0.5 ; sigmaWTrue <- 1/6
sigmaEpsTrue <- 0.35
sigmaXTrue <- 0.1

# Simulate data with missingess not at random
# according to a linear logistic regression model:

set.seed(2)
n <- 1000
w <- rnorm(n,muWTrue,sigmaWTrue)
x <- rnorm(n,w,sigmaXTrue)
y <- fTrue(x) + sigmaEpsTrue*rnorm(n)

# Split into observed and unobserved data according
# to specified proportion:

propObs <- 0.10
nObs <- round(propObs*n)
nUnobs <- n - nObs
indsObs <- sample(1:n,nObs,replace = FALSE)
indsUnobs <- setdiff(1:n,indsObs)

# Obtain observed data:

xObs <-  x[indsObs]
wxObs <- w[indsObs]
yxObs <- y[indsObs]

# Obtain unobserved data:

xUnobs <- x[indsUnobs]
wxUnobs <- w[indsUnobs]
yxUnobs <- y[indsUnobs]

# Obtain the Z matrix for the observed data:

ncZ <- 35
knots <- seq(min(xObs),max(xObs),length = (ncZ+2))[-c(1,(ncZ+2))]
ZxObs <- outer(xObs,knots,"-")
ZxObs <- ZxObs*(ZxObs>0)

# Specify model in Stan:

npRegBerksonMeaErrModel <- 
   'data
   {
      int<lower=1> nObs;        int<lower=1> nUnobs;
      int<lower=1> ncZ;
      vector[nObs] yxObs;       vector[nUnobs] yxUnobs;
      vector[nObs] wxObs;       vector[nUnobs] wxUnobs;
      vector[nObs] xObs;        
      vector[ncZ]  knots;       matrix[nObs,ncZ] ZxObs;
      real<lower=0> sigmaBeta;  real<lower=0> sigmaMu;
      real<lower=0> Ax;         real<lower=0> Aeps;
      real<lower=0> Au;           
   }
   parameters 
   {
      vector[2] beta;           vector[ncZ] u;
      real<lower=0> sigmaX;     real<lower=0> sigmaEps;
      real<lower=0> sigmaU;     real xUnobs[nUnobs];
   }
   transformed parameters 
   {
      matrix[nObs,2] XxObs;     matrix[nUnobs,2] XxUnobs;
      matrix[nUnobs,ncZ] ZxUnobs;
      for (i in 1:nObs)
      {
         XxObs[i,1] = 1    ;   XxObs[i,2] = xObs[i];
      }
      for (i in 1:nUnobs)
      {
         XxUnobs[i,1] = 1    ;   XxUnobs[i,2] = xUnobs[i];
         for (k in 1:ncZ)   
            ZxUnobs[i,k] = (xUnobs[i]-knots[k])*step(xUnobs[i]-knots[k]);
      }
   }
   model 
   {
      yxObs ~ normal(XxObs*beta+ZxObs*u,sigmaEps);
      xObs  ~ normal(wxObs,sigmaX); 

      yxUnobs ~ normal(XxUnobs*beta+ZxUnobs*u,sigmaEps); 
      xUnobs  ~ normal(wxUnobs,sigmaX); 

      u ~ normal(0,sigmaU) ;   beta ~ normal(0,sigmaBeta);
      sigmaX ~ cauchy(0,Ax);   sigmaEps ~ cauchy(0,Aeps);   
      sigmaU ~ cauchy(0,Au);
   }'

# Set up input data:

allData <- list(nObs = nObs,nUnobs = nUnobs,ncZ = ncZ,xObs = xObs,wxObs = wxObs,
                wxUnobs = wxUnobs,yxObs = yxObs,yxUnobs = yxUnobs,ZxObs = ZxObs,
                knots = knots,sigmaMu = 1e5,sigmaBeta = 1e5,
                Ax = 1e5,Aeps = 1e5,Au = 1e5)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = npRegBerksonMeaErrModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = npRegBerksonMeaErrModel,data = allData,warmup = nWarm,
                 iter = (nWarm  + nKept),chains = 1,thin = nThin,refresh = 25,
                 fit = stanCompilObj)

# Extract relevant MCMC samples:

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

xyObsCol <- "darkblue" 
wyCol <- "lightskyblue"
cexVal <- 0.3
xlimVal <- c(0.1,0.9)

estFunCol <- "darkgreen"; trueFunCol <- "indianred3"
varBandCol <- "palegreen"
ng <- 201
xg <- seq(xlimVal[1],xlimVal[2],length = ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- outer(xg,knots,"-")
Zg <- Zg*(Zg>0)
fMCMC <-   Xg%*%betaMCMC + Zg%*%uMCMC 
credLower <- apply(fMCMC,1,quantile,0.025)
credUpper <- apply(fMCMC,1,quantile,0.975)
fg <- apply(fMCMC,1,mean)

par(mfrow = c(1,1),mai = c(1.02,0.82,0.82,0.42))

plot(xg,fg,type = "n",bty = "l",xlab = "x",ylab = "y",xlim = xlimVal,ylim = range(y),
     cex.lab = 1.5)
polygon(c(xg,rev(xg)),c(credLower,rev(credUpper)),col = varBandCol,border = FALSE)

lines(xg,fg,lwd = 2,col = estFunCol)
lines(xg,fTrue(xg),lwd = 2,col = trueFunCol)
points(xObs,yxObs,col = xyObsCol,cex = cexVal)
points(wxObs,yxObs,col = wyCol,cex = cexVal)
points(wxUnobs,yxUnobs,col = wyCol,cex = cexVal)
legend("bottomleft",c("validation data",
       "data with predictor subject to measurement error"),
       col = c(xyObsCol,wyCol),pch = rep(1,2),pt.cex = cexVal/1.3,
       cex = 1.3)

legend("topright",c("true f","estimated f"),col = c(trueFunCol,estFunCol),
       lwd = rep(2,2),cex = 1.3)

# Display summary of MCMC samples for model parameters:

indQ1 <- length(xg[xg <= quantile(xObs,0.25)])
fQ1MCMC <- fMCMC[indQ1,]
indQ2 <- length(xg[xg <= quantile(xObs,0.50)])
fQ2MCMC <- fMCMC[indQ2,]
indQ3 <- length(xg[xg <= quantile(xObs,0.75)])
fQ3MCMC <- fMCMC[indQ3,]

parms <- list(cbind(sigmaEpsMCMC,sigmaXMCMC,fQ1MCMC,fQ2MCMC,fQ3MCMC))
parNamesVal <- list(c(expression(sigma[epsilon])),c(expression(sigma[x])),
                      c("funct. est. at 1st","quantile of",
                        "obs'd predictor"),
                      c("funct. est. at 2nd","quantile of",
                        "obs'd predictor"),
                      c("funct. est. at 3rd","quantile of",
                        "obs'd predictor"))

summMCMC(parms,parNames = parNamesVal,numerSummCex = 1.25,
         columnHeadCex = 2.9,KDEvertLine = FALSE,
         addTruthToKDE = c(sigmaEpsTrue,sigmaXTrue,fTrue(xg[indQ1]),
         fTrue(xg[indQ2]),fTrue(xg[indQ3])))

# Display summary for some of the missing x values:

parms <- list(t(xUnobsMCMC[1:5,]))
parNamesVal <- list(c(expression(x["unobs,1"])),
                      c(expression(x["unobs,2"])),
                      c(expression(x["unobs,3"])),
                      c(expression(x["unobs,4"])),
                      c(expression(x["unobs,5"])))

summMCMC(parms,parNames = parNamesVal,KDEvertLine = FALSE,,columnHeadCex = 2.9,
         addTruthToKDE = x[1:5])

########## End of npRegBerksonMeaErr ##########
