########## R script: indonRespirBayes ##########

# For conducting a Bayesian additive mixed
# model analysis on the spinal bone mineral
# density data using Markov chain Monte Carlo 
# and Stan.

# Last changed: 22 AUG 2017

# Set flag for code compilation (needed if 
# running script first time in current session):

compileCode <- TRUE

# Set MCMC sample size parameters:

nWarm <- 5000
nKept <- 5000
nThin <- 5

# Load required packages:

library(HRW)   ;   library(rstan)   ;   library(lattice)

# Load in the Indonesian respiratory data and extract the 
# component variables:

data(indonRespir)

idnum <- indonRespir$idnum
y <- indonRespir$respirInfec
x1 <- indonRespir$vitAdefic
x2 <- 1 - indonRespir$female
x3Orig <- indonRespir$height
x4 <- indonRespir$stunted
x5 <- indonRespir$visit2
x6 <- indonRespir$visit3
x7 <- indonRespir$visit4
x8 <- indonRespir$visit5
x9 <- indonRespir$visit6
x10Orig <- indonRespir$age

numObs <- length(y)
numGrp <- length(unique(idnum))

# Standardize continuous data for Bayesian analysis:

mean.x3 <- mean(x3Orig)     ;   sd.x3 <- sd(x3Orig)
mean.x10 <- mean(x10Orig)   ;   sd.x10 <- sd(x10Orig)
x3 <- (x3Orig - mean.x3)/sd.x3
x10 <- (x10Orig - mean.x10)/sd.x10

# Set up matrices for additive model:

numObs <- length(y)
X <- cbind(rep(1,numObs),x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
ncX <- ncol(X)

numIntKnots <- 15
intKnots <- quantile(unique(x10),seq(0,1,length = 
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])
range.x10 <- c(1.01*min(x10)-0.01*max(x10),1.01*max(x10)-0.01*min(x10))
Zspl <- ZOSull(x10,intKnots = intKnots,range.x = range.x10)
ncZspl <- ncol(Zspl)

numSpl <- ncol(Zspl)
numGrp <- length(unique(idnum))
numObs <- length(y)

# Set hyperparameters:

sigmaBeta <- 1e5     ;     Agrp <- 1e5    ;     Aspl <- 1e5

# Specify model in Stan:

logistAddMixModModel <- 
'data
{
   int<lower=1> numObs;              int<lower=1> numGrp;
   int<lower=1> ncX;                 int<lower=1> ncZspl;        
   real<lower=0> sigmaBeta;
   real<lower=0> Agrp;               real<lower=0> Aspl;
   int<lower=0,upper=1> y[numObs];   int<lower=1> idnum[numObs];
   matrix[numObs,ncX] X;             matrix[numObs,ncZspl] Zspl;
}
parameters
{
   vector[ncX] beta;            vector[numGrp] U;            
   vector[ncZspl] u;            real<lower=0> sigmaGrp;      
   real<lower=0> sigmaSpl;  
}
model
{
   y ~ bernoulli_logit(X*beta + U[idnum] + Zspl*u);	
   U ~ normal(0,sigmaGrp);           u ~ normal(0,sigmaSpl);
   beta  ~ normal(0,sigmaBeta) ;     sigmaGrp ~ cauchy(0,Agrp) ;
   sigmaSpl ~ cauchy(0,Aspl);
}'

# Fit model using MCMC via Stan:

allData <- list(numObs = numObs,numGrp = numGrp,ncX = ncX,ncZspl = ncZspl,
                idnum = idnum,X = X,y = y,Zspl = Zspl,sigmaBeta = sigmaBeta,
                Agrp = Agrp,Aspl = Aspl)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = logistAddMixModModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = logistAddMixModModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 1000,
                 fit = stanCompilObj)

# Save and extract relevant MCMC samples:

beta0MCMC <- extract(stanObj,"beta[1]",permuted = FALSE)
beta1MCMC <- extract(stanObj,"beta[2]",permuted = FALSE)
beta2MCMC <- extract(stanObj,"beta[3]",permuted = FALSE)
beta3MCMC <- extract(stanObj,"beta[4]",permuted = FALSE)
beta4MCMC <- extract(stanObj,"beta[5]",permuted = FALSE)
beta5MCMC <- extract(stanObj,"beta[6]",permuted = FALSE)
beta6MCMC <- extract(stanObj,"beta[7]",permuted = FALSE)
beta7MCMC <- extract(stanObj,"beta[8]",permuted = FALSE)
beta8MCMC <- extract(stanObj,"beta[9]",permuted = FALSE)
beta9MCMC <- extract(stanObj,"beta[10]",permuted = FALSE)
beta10MCMC <- extract(stanObj,"beta[11]",permuted = FALSE)
sigmaGrpMCMC <- extract(stanObj,"sigmaGrp",permuted = FALSE)
sigmaSplMCMC <- extract(stanObj,"sigmaSpl",permuted = FALSE)

uMCMC <- NULL
for (k in 1:ncZspl)
{
   charVar <- paste("u[",as.character(k),"]",sep = "") 
   uMCMC <- rbind(uMCMC,extract(stanObj,charVar,permuted = FALSE))
}

# Convert to coefficient for continuous predictor to correspond
# to original scale:

beta3MCMCorig <- beta3MCMC/sd.x3

# Do parameters plot:

parms <- list(cbind(beta1MCMC,beta2MCMC,beta3MCMCorig,beta4MCMC,beta5MCMC,
                    beta6MCMC,beta7MCMC,beta8MCMC,beta9MCMC,sigmaGrpMCMC))
parNamesVal <-   list(c("vitamin A","deficiency"),c("male"),c("height"),
                      c("stunted"), c("2 visits"),c("3 visits"),
                      c("4 visits"),c("5 visits"),c("6 visits"),
                      c(expression(sigma[grp])))
  
summMCMC(parms,parNames = parNamesVal)

# Do plot for estimated probability function given age:

par(mai = c(1,1.2,0.5,0.2))
cex.labVal <- 1.8   ;  cex.axisVal <- 1.5
ng <- 101
x10g <- seq(min(x10),max(x10),length = ng)
XotherPreds <- X[,2:10]
otherPredsMeans <- apply(XotherPreds,2,mean)
covarAddOn <- t(matrix(rep(otherPredsMeans,ng),
                    length(otherPredsMeans),ng))
Xg <- cbind(rep(1,ng),covarAddOn,x10g)
Zg <- ZOSull(x10g,intKnots = intKnots,range.x = range.x10)

betaMCMC <- rbind(beta0MCMC,beta1MCMC,beta2MCMC,beta3MCMCorig,beta4MCMC,beta5MCMC,
                  beta6MCMC,beta7MCMC,beta8MCMC,beta9MCMC,beta10MCMC)

etaHatMCMC <- Xg%*%betaMCMC + Zg%*%uMCMC
muHatMCMC <- 1/(1+exp(-etaHatMCMC))

credLower <- apply(muHatMCMC,1,quantile,0.025)
credUpper <- apply(muHatMCMC,1,quantile,0.975)
muHatg <- apply(muHatMCMC,1,mean)
ylimVal <- range(c(credLower,credUpper))

x10gOrig <- sd.x10*x10g + mean.x10
plot(x10gOrig,muHatg,type = "n",xlab = "age in years",
     ylab = "estimated probability of respiratory infection",
     ylim = ylimVal,bty = "l",cex.lab = cex.labVal,cex.axis = cex.axisVal)
polygon(c(x10gOrig,rev(x10gOrig)),c(credLower,rev(credUpper)),
        col = "palegreen",border = FALSE)
lines(x10gOrig,muHatg,col = "darkgreen",lwd = 2)
rug(jitter(x10Orig),col = "dodgerblue")

# Do a summaries of the MCMC samples for slices of
# the logit probability function at the quartiles:

indQ1 <- length(x10g[x10g<quantile(x10,0.25)])
indQ2 <- length(x10g[x10g<quantile(x10,0.50)])
indQ3 <- length(x10g[x10g<quantile(x10,0.75)])

etaHatMCMCQ1 <- etaHatMCMC[indQ1,]
etaHatMCMCQ2 <- etaHatMCMC[indQ2,]
etaHatMCMCQ3 <- etaHatMCMC[indQ3,]

parms <- list(cbind(etaHatMCMCQ1,etaHatMCMCQ2,etaHatMCMCQ3))
parNamesVal <- list(c("logit probab.","respir. infec.","at 1st quart. age"),
                    c("logit probab.","respir. infec.","at 2nd quart. age"),
                    c("logit probab.","respir. infec.","at 3rd quart. age"))
summMCMC(parms,parNames = parNamesVal)

########## End of indonRespirBayes ##########


