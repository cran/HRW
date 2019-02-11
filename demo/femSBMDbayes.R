########## R script: femSBMDbayes ##########

# For conducting a Bayesian additive mixed
# model analysis on the spinal bone mineral
# density data using Markov chain Monte Carlo 
# and Stan.

# Last changed: 29 AUG 2018

# Set flag for code compilation (needed if 
# running script first time in current session):

compileCode <- TRUE

# Set MCMC sample size paramaters:

nWarm <- 1000
nKept <- 1000
nThin <- 1

# Load required packages:

library(HRW)   ;   library(rstan)   ;   library(lattice)

# Load in the female spinal bone data and extract the 
# component variables:

data(femSBMD)
idnum <- femSBMD$idnum 
x1 <- femSBMD$black 
x2 <- femSBMD$hispanic 
x3 <- femSBMD$white
x4Orig <- femSBMD$age 
yOrig <- femSBMD$spnbmd

# Standardise data for Bayesian analysis:

mean.x4 <- mean(x4Orig) ; sd.x4 <- sd(x4Orig)
mean.y <- mean(yOrig)   ; sd.y <- sd(yOrig)
x4 <- (x4Orig - mean.x4)/sd.x4
y <- (yOrig - mean.y)/sd.y

# Set up matrices for additive model:

numObs <- length(y)
X <- cbind(rep(1,numObs),x1,x2,x3,x4)
ncX <- ncol(X)

numIntKnots <- 15
intKnots <- quantile(unique(x4),seq(0,1,length = 
                 (numIntKnots+2))[-c(1,(numIntKnots+2))])
range.x4 <- c(1.01*min(x4)-0.01*max(x4),1.01*max(x4)-0.01*min(x4))
Zspl <- ZOSull(x4,intKnots = intKnots,range.x = range.x4)
ncZ <- ncol(Zspl)

numGrp <- length(unique(idnum))
numObs <- length(y)

# Set hyperparameters:

sigmaBeta <- 1e5 ; AU <- 1e5 ; Au <- 1e5 ; Aeps <- 1e5

# Specify model in Stan:

addMixModModel <- 
'data
{
   int<lower=1> numObs;         int<lower=1> numGrp;
   int<lower=1> ncX;            int<lower=1> ncZ;        
   real<lower=0> sigmaBeta;     real<lower=0> AU;
   real<lower=0> Aeps;          real<lower=0> Au;
   vector[numObs] y;            int<lower=1> idnum[numObs];
   matrix[numObs,ncX] X;        matrix[numObs,ncZ] Zspl;
}
parameters
{
   vector[ncX] beta;            vector[numGrp] U;            
   vector[ncZ] u;               real<lower=0> sigmaU;        
   real<lower=0> sigmau;        real<lower=0> sigmaEps;     
}
model
{
   y ~ normal(X*beta + U[idnum] + Zspl*u,sigmaEps);	
   U ~ normal(0,sigmaU);           u ~ normal(0,sigmau);
   beta  ~ normal(0,sigmaBeta) ;   sigmaEps ~ cauchy(0,Aeps);
   sigmaU ~ cauchy(0,AU) ;   sigmau ~ cauchy(0,Au);
}'

# Fit model using MCMC via Stan:

allData <- list(numObs = numObs,numGrp = numGrp,ncX = ncX,ncZ = ncZ,
                idnum = idnum,X = X,y = y,Zspl = Zspl,sigmaBeta = sigmaBeta,
                AU = AU,Au = Au,Aeps = Aeps)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = addMixModModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = addMixModModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 100,
                 fit = stanCompilObj)

# Save and extract relevant MCMC samples:

beta0MCMC <- as.vector(extract(stanObj,"beta[1]",permuted = FALSE))
beta1MCMC <- as.vector(extract(stanObj,"beta[2]",permuted = FALSE))
beta2MCMC <- as.vector(extract(stanObj,"beta[3]",permuted = FALSE))
beta3MCMC <- as.vector(extract(stanObj,"beta[4]",permuted = FALSE))
beta4MCMC <- as.vector(extract(stanObj,"beta[5]",permuted = FALSE))

sigmaUMCMC <- as.vector(extract(stanObj,"sigmaU",permuted = FALSE))
sigmaEpsMCMC <- as.vector(extract(stanObj,"sigmaEps",permuted = FALSE))


UMCMC <- NULL
for (k in 1:numGrp)
{
   charVar <- paste("U[",as.character(k),"]",sep = "") 
   UMCMC <- rbind(UMCMC,extract(stanObj,charVar,permuted = FALSE))
}

uMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("u[",as.character(k),"]",sep = "") 
   uMCMC <- rbind(uMCMC,extract(stanObj,charVar,permuted = FALSE))
}

# Convert to parameters original scale:

beta1MCMCorig <- beta1MCMC*sd.y
beta2MCMCorig <- beta2MCMC*sd.y
beta3MCMCorig <- beta3MCMC*sd.y
beta4MCMCorig <- beta4MCMC*(sd.y/sd.x4)

beta0MCMCorig <-  mean.y + sd.y*beta0MCMC -  mean.x4*beta4MCMCorig

sigmaUMCMCorig <- sigmaUMCMC*sd.y
sigmaEpsMCMCorig <- sigmaEpsMCMC*sd.y

# Do parameters plot:

parms <- list(cbind(beta1MCMCorig,beta2MCMCorig,beta3MCMCorig,
                    sigmaUMCMCorig,sigmaEpsMCMCorig))
parNamesVal <- list(c("Black","vs. Asian"),c("Hispanic","vs. Asian"),
                      c("White","vs. Asian"),c(expression(sigma[U])),
                      c(expression(sigma[epsilon])))
summMCMC(parms,parNames = parNamesVal)

# Do fitted curves using lattice graphics:

ng <- 101
ylim.val <- c(min(yOrig),max(yOrig))
x4g <- seq(min(x4),max(x4),length = ng)
Xg <- cbind(rep(1,ng),x4g)
Zg <- ZOSull(x4g,intKnots = intKnots,range.x = range.x4)

midCurvsg <- list()
lowCurvsg <- list()
uppCurvsg <- list()

for (ip in 1:4)
{
   if (ip == 1) betaSplMCMC <- rbind(beta0MCMC,beta4MCMC)
   if (ip == 2) betaSplMCMC <- rbind(beta0MCMC + beta1MCMC,beta4MCMC)
   if (ip == 3) betaSplMCMC <- rbind(beta0MCMC + beta2MCMC,beta4MCMC)
   if (ip == 4) betaSplMCMC <- rbind(beta0MCMC + beta3MCMC,beta4MCMC)

   fMCMC <-   Xg%*%betaSplMCMC + Zg%*%uMCMC
   fMCMCorig <- fMCMC*sd.y + mean.y
   credLowerOrig <- apply(fMCMCorig,1,quantile,0.025)
   credUpperOrig <- apply(fMCMCorig,1,quantile,0.975)
   fhatgOrig <- apply(fMCMCorig,1,mean)

   x4gOrig <- x4g*sd.x4 + mean.x4

   midCurvsg[[ip]] <- apply(fMCMCorig,1,mean)
   lowCurvsg[[ip]] <- apply(fMCMCorig,1,quantile,0.025)
   uppCurvsg[[ip]] <- apply(fMCMCorig,1,quantile,0.975)
}

fitFig <- xyplot(yOrig~x4Orig|ethnicity,groups = idnum,
                  as.table = TRUE,data = femSBMD,
                  strip = strip.custom(par.strip.text = list(cex = 1.5)),
                  par.settings = list(layout.heights = list(strip = 1.6)),
                  scales = list(cex = 1.25),
                  xlab = list("age (years)",cex = 1.5),
                  ylab = list(expression(paste(
                  "spinal bone mineral density (g/c",m^2,")")),cex = 1.5),
                  subscripts = TRUE,
                  panel = function(x,y,subscripts,groups)
                  {
                     panel.grid() 
                     if (any(femSBMD$ethnicity[subscripts] == "Asian"))    panNum <- 1
                     if (any(femSBMD$ethnicity[subscripts] == "Black"))    panNum <- 2
                     if (any(femSBMD$ethnicity[subscripts] == "Hispanic")) panNum <- 3
                     if (any(femSBMD$ethnicity[subscripts] == "White"))    panNum <- 4

                     panel.superpose(x,y,subscripts,groups,type = "b",pch = 16)

                     panel.polygon(c(x4gOrig,rev(x4gOrig)),
                                 c(lowCurvsg[[panNum]],rev(uppCurvsg[[panNum]])),
                                 col = "palegreen",border = FALSE)

                     panel.xyplot(x4gOrig,midCurvsg[[panNum]],lwd = 2,type = "l",col = "darkgreen")     
                })
   
print(fitFig)

########## End of femSBMDbayes ##########


