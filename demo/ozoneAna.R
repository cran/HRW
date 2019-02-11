######### R script: ozoneAna ##########

# For fitting a penalized spline-based additive 
# to the ozone data when the predictors are subject 
# to missingness completely at random.

# Last changed: 22 AUG 2017

# Set flag for code compilation (needed if 
# running script first time in current session):

compileCode <- TRUE

# Load required packages:

library(rstan)
library(HRW)
library(mlbench)

# Set MCMC parameters:

nWarm <- 1000          # Length of burn-in.
nKept <- 1000            # Size of the kept sample.
nThin <- 1               # Thinning factor. 

# Obtain data:

data(Ozone)
names(Ozone) <- c("month","day","week","maxOzone","pressHt",
                  "windSpeed","humidity","tempSandberg",
                  "tempElMonte","inverseBaseHt","daggettPresGrad",
                  "inverseBaseTemp","visibility")

# Set up predictors:

rx1 <- as.numeric(!is.na(Ozone$tempElMonte))
x1ObsOrig <- Ozone$tempElMonte[rx1 == 1]
mean.x1Obs <- mean(x1ObsOrig)  ; sd.x1Obs <- sd(x1ObsOrig)
x1Obs <- (x1ObsOrig - mean.x1Obs)/sd.x1Obs

rx2 <- as.numeric(!is.na(Ozone$daggettPresGrad))
x2ObsOrig <- Ozone$daggettPresGrad[rx2 == 1]
mean.x2Obs <- mean(x2ObsOrig)  ; sd.x2Obs <- sd(x2ObsOrig)
x2Obs <- (x2ObsOrig - mean.x2Obs)/sd.x2Obs

rx3 <- as.numeric(!is.na(Ozone$inverseBaseHt))
x3ObsOrig <- Ozone$inverseBaseHt[rx3 == 1]
mean.x3Obs <- mean(x3ObsOrig)  ; sd.x3Obs <- sd(x3ObsOrig)
x3Obs <- (x3ObsOrig - mean.x3Obs)/sd.x3Obs

# Set up response 

ry <- as.numeric(!is.na(Ozone$maxOzone))
yObsOrig <- Ozone$maxOzone[ry == 1]
mean.yObs <- mean(yObsOrig)   ;   sd.yObs <- sd(yObsOrig)
yObs <- (yObsOrig - mean.yObs)/sd.yObs

# Set sample sizes:

n <- nrow(Ozone)
nObsy <- length(yObs)     ;   nObsx1 <- length(x1Obs)  
nObsx2 <- length(x2Obs)   ;   nObsx3 <- length(x3Obs)  

# Obtain the knots based on the observed data:

ncZ1 <- 30
knots1 <- seq(min(x1Obs),max(x1Obs),length = (ncZ1+2))[-c(1,ncZ1+2)]
Z1Obs <- outer(x1Obs,knots1,"-")
Z1Obs <- Z1Obs*(Z1Obs>0)

ncZ2 <- 30
knots2 <- seq(min(x2Obs),max(x2Obs),length = (ncZ2+2))[-c(1,ncZ2+2)]
Z2Obs <- outer(x2Obs,knots2,"-")
Z2Obs <- Z2Obs*(Z2Obs>0)

ncZ3 <- 30
knots3 <- seq(min(x3Obs),max(x3Obs),length = (ncZ3+2))[-c(1,ncZ3+2)]
Z3Obs <- outer(x3Obs,knots3,"-")
Z3Obs <- Z3Obs*(Z3Obs>0)

# Specify model in Stan:

ozoneAddModMCARModel <- 
   'data
   {
      int<lower=1> n;
      int<lower=1> nObsy;         int<lower=1> nObsx1; 
      int<lower=1> nObsx2;        int<lower=1> nObsx3;        
      int<lower=1> ncZ1;          int<lower=1> ncZ2;         
      int<lower=1> ncZ3;           
      vector[nObsy] yObs;        vector[nObsx1] x1Obs;       
      vector[nObsx2] x2Obs;      vector[nObsx3] x3Obs;        
      vector[ncZ1]  knots1;      vector[ncZ2]  knots2; 
      vector[ncZ3]  knots3;       
      matrix[nObsx1,ncZ1] Z1Obs;  matrix[nObsx2,ncZ2] Z2Obs;
      matrix[nObsx3,ncZ3] Z3Obs;  
      real<lower=0> sigmaBeta;   real<lower=0> sigmaGamma;
      real<lower=0> sigmaMu;     real<lower=0> Ax;         
      real<lower=0> Au;          real<lower=0> Av;
   }
   parameters 
   {
      vector[4] beta;            vector[4] gamma;
      vector[ncZ1] u1;           vector[ncZ1] v1;
      vector[ncZ2] u2;           vector[ncZ2] v2;
      vector[ncZ3] u3;           vector[ncZ3] v3;
      vector[3] muX;             cov_matrix[3] SigmaX;
      real<lower=0> sigmaU1;     real<lower=0> sigmaV1;
      real<lower=0> sigmaU2;     real<lower=0> sigmaV2;
      real<lower=0> sigmaU3;     real<lower=0> sigmaV3;
      vector[n-nObsy] yUnobs;    vector[n-nObsx1] x1Unobs;  
      vector[n-nObsx2] x2Unobs;  vector[n-nObsx3] x3Unobs;  
      vector[3] a; 
   }
   transformed parameters 
   {
      vector[n] y;        vector[n] x1;        
      vector[n] x2;       vector[n] x3;        
      matrix[3,n] xMat;   matrix[n,ncZ1] Z1;
      matrix[n,ncZ2] Z2;  matrix[n,ncZ3] Z3;  
      vector[n] meanVec;     vector[n] logvar;

      for (i in 1:nObsy) y[i] = yObs[i];
      for (i in 1:(n-nObsy)) y[i+nObsy] = yUnobs[i];

      for (i in 1:nObsx1) x1[i] = x1Obs[i];
      for (i in 1:(n-nObsx1)) x1[i+nObsx1] = x1Unobs[i];

      for (i in 1:nObsx2) x2[i] = x2Obs[i];
      for (i in 1:(n-nObsx2)) x2[i+nObsx2] = x2Unobs[i];

      for (i in 1:nObsx3) x3[i] = x3Obs[i];
      for (i in 1:(n-nObsx3)) x3[i+nObsx3] = x3Unobs[i];

      for (i in 1:n) 
      {
         xMat[1,i] = x1[i] ;  xMat[2,i] = x2[i] ; 
         xMat[3,i] = x3[i] ;   
      }

      for (k in 1:ncZ1)
      {           
         for (i in 1:nObsx1) Z1[i,k] = Z1Obs[i,k];
         for (i in (nObsx1+1):n) 
            Z1[i,k] = (x1[i]-knots1[k])*step(x1[i]-knots1[k]);
      }
      for (k in 1:ncZ2)
      {           
         for (i in 1:nObsx2) Z2[i,k] = Z2Obs[i,k];
         for (i in (nObsx2+1):n) 
            Z2[i,k] = (x2[i]-knots2[k])*step(x2[i]-knots2[k]);
      }
      for (k in 1:ncZ3)
      {           
         for (i in 1:nObsx3) Z3[i,k] = Z3Obs[i,k];
         for (i in (nObsx3+1):n) 
            Z3[i,k] = (x3[i]-knots3[k])*step(x3[i]-knots3[k]);
      }
 
      for (i in 1:n)
      {
         meanVec[i] = beta[1] + beta[2]*x1[i] + beta[3]*x2[i] 
                            + beta[4]*x3[i] 
                            + dot_product(u1,Z1[i]) 
                            + dot_product(u2,Z2[i]) 
                            + dot_product(u3,Z3[i]) ;

         logvar[i]= gamma[1] + gamma[2]*x1[i] + gamma[3]*x2[i] 
                              + gamma[4]*x3[i] 
                              + dot_product(v1,Z1[i]) 
                              + dot_product(v2,Z2[i]) 
                              + dot_product(v3,Z3[i]) ;
      }
   }
   model 
   {
      vector[3] zeroVec;
      matrix[3,3] scaleMatrix;

      y ~ normal(meanVec,exp(logvar)); 
      u1 ~ normal(0,sigmaU1) ;      v1 ~ normal(0,sigmaV1) ;
      u2 ~ normal(0,sigmaU2) ;      v2 ~ normal(0,sigmaV2) ;
      u3 ~ normal(0,sigmaU3) ;      v3 ~ normal(0,sigmaV3) ;
      beta ~ normal(0,sigmaBeta);   gamma ~ normal(0,sigmaGamma);
      sigmaU1 ~ cauchy(0,Au);       sigmaV1 ~ cauchy(0,Av);
      sigmaU2 ~ cauchy(0,Au);       sigmaV2 ~ cauchy(0,Av);
      sigmaU3 ~ cauchy(0,Au);       sigmaV3 ~ cauchy(0,Av);

      for (i in 1:n)
        col(xMat,i) ~ multi_normal(muX,SigmaX);

      for (j in 1:3)
      {
         muX[j] ~ normal(0,Ax);
         for (jd in (j+1):3)
         {
            scaleMatrix[j,jd] = 0;     
            scaleMatrix[jd,j] = 0;     
         }
         a[j] ~ inv_gamma(0.5,pow(Ax,-2));     
         scaleMatrix[j,j] = 4/a[j];
      }  
      SigmaX ~ inv_wishart(4,scaleMatrix);
   }'

# Set up input data:

allData <- list(n=n,yObs=yObs,nObsy=nObsy,
                nObsx1 = nObsx1,nObsx2 = nObsx2,nObsx3 = nObsx3,
                x1Obs = x1Obs,x2Obs = x2Obs,x3Obs = x3Obs,
                ncZ1 = ncZ1,ncZ2 = ncZ2,ncZ3 = ncZ3,
                knots1 = knots1,knots2 = knots2,knots3 = knots3,
                Z1Obs = Z1Obs,Z2Obs = Z2Obs,Z3Obs = Z3Obs,
                sigmaMu = 1e5,sigmaBeta = 1e5,sigmaGamma = 1e5,
                Ax = 1e5,Au = 1e5,Av = 1e5)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = ozoneAddModMCARModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = ozoneAddModMCARModel,data = allData,
                 warmup = nWarm,iter = (nWarm + nKept),chains = 1,
                 thin = nThin,refresh = 10,fit = stanCompilObj)

# Extract relevant MCMC sampes:

betaMCMC <- NULL  ; gammaMCMC <- NULL
for (j in 1:4)
{
   charVar <- paste("beta[",as.character(j),"]",sep = "") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar,permuted = FALSE))
   charVar <- paste("gamma[",as.character(j),"]",sep = "") 
   gammaMCMC <- rbind(gammaMCMC,extract(stanObj,charVar,permuted = FALSE))
}

u1MCMC <- NULL  ;  v1MCMC <- NULL
for (k in 1:ncZ1)
{
   charVar <- paste("u1[",as.character(k),"]",sep = "") 
   u1MCMC <- rbind(u1MCMC,extract(stanObj,charVar,permuted = FALSE))
   charVar <- paste("v1[",as.character(k),"]",sep = "") 
   v1MCMC <- rbind(v1MCMC,extract(stanObj,charVar,permuted = FALSE))
}
u2MCMC <- NULL  ;  v2MCMC <- NULL
for (k in 1:ncZ2)
{
   charVar <- paste("u2[",as.character(k),"]",sep = "") 
   u2MCMC <- rbind(u2MCMC,extract(stanObj,charVar,permuted = FALSE))
   charVar <- paste("v2[",as.character(k),"]",sep = "") 
   v2MCMC <- rbind(v2MCMC,extract(stanObj,charVar,permuted = FALSE))
}
u3MCMC <- NULL  ;  v3MCMC <- NULL
for (k in 1:ncZ3)
{
   charVar <- paste("u3[",as.character(k),"]",sep = "") 
   u3MCMC <- rbind(u3MCMC,extract(stanObj,charVar,permuted = FALSE))
   charVar <- paste("v3[",as.character(k),"]",sep = "") 
   v3MCMC <- rbind(v3MCMC,extract(stanObj,charVar,permuted = FALSE))
}
uMCMC <- rbind(u1MCMC,u2MCMC,u3MCMC)
vMCMC <- rbind(v1MCMC,v2MCMC,v3MCMC)

yUnobsMCMC <- NULL
for (i in 1:(n-nObsy))
{
   charVar <- paste("yUnobs[",as.character(i),"]",sep = "") 
   yUnobsMCMC <- rbind(yUnobsMCMC,extract(stanObj,charVar,permuted = FALSE))
}

x1UnobsMCMC <- NULL
for (i in 1:(n-nObsx1))
{
   charVar <- paste("x1Unobs[",as.character(i),"]",sep = "") 
   x1UnobsMCMC <- rbind(x1UnobsMCMC,extract(stanObj,charVar,permuted = FALSE))
}

x2UnobsMCMC <- NULL
for (i in 1:(n-nObsx2))
{
   charVar <- paste("x2Unobs[",as.character(i),"]",sep = "") 
   x2UnobsMCMC <- rbind(x2UnobsMCMC,extract(stanObj,charVar,permuted = FALSE))
}

x3UnobsMCMC <- NULL
for (i in 1:(n-nObsx3))
{
   charVar <- paste("x3Unobs[",as.character(i),"]",sep = "") 
   x3UnobsMCMC <- rbind(x3UnobsMCMC,extract(stanObj,charVar,permuted = FALSE))
}

# Plot fit:

par(mfcol = c(2,3),mai = c(0.8,0.7,0.08,0.03))

# Plot first mean function fit (with other variables
# set at their averages):

cex.labVal <- 1.7
ng <- 101 ; 
ylimValMean <- c(2,20) ; ylimValStDev <- c(2,11)
x1g <- seq(min(x1Obs),max(x1Obs),length = ng)
x2g <- rep(mean(x2Obs),length = ng)
x3g <- rep(mean(x3Obs),length = ng)
Z1g <- outer(x1g,knots1,"-") ; Z1g <- Z1g*(Z1g>0)
Z2g <- outer(x2g,knots2,"-") ; Z2g <- Z2g*(Z2g>0)
Z3g <- outer(x3g,knots3,"-") ; Z3g <- Z3g*(Z3g>0)
Xg <- cbind(rep(1,ng),x1g,x2g,x3g)
Zg <- cbind(Z1g,Z2g,Z3g)
f1MCMC <- Xg%*%betaMCMC + Zg%*%uMCMC
f1MCMCOrig <- sd.yObs*f1MCMC + mean.yObs
credLower <- apply(f1MCMCOrig,1,quantile,0.025)
credUpper <- apply(f1MCMCOrig,1,quantile,0.975)
f1gOrig <- apply(f1MCMCOrig,1,mean)
x1gOrig <- sd.x1Obs*x1g + mean.x1Obs
plot(x1gOrig,f1gOrig,type = "n",bty = "l",
     xlab = expression(paste("El Monte temperature (",phantom()^{o},"F)")),
     ylab = "max. 1-hour-average ozone level",ylim = ylimValMean,
     cex.lab = cex.labVal)
polygon(c(x1gOrig,rev(x1gOrig)),c(credLower,rev(credUpper)),
        col = "palegreen",border = FALSE)
lines(x1gOrig,f1gOrig,col = "darkgreen",lwd = 2)
rug(x1ObsOrig,col = "dodgerblue")
legend("topleft","mean functions",cex = 1.4)

# Plot first standard deviation function fit (with other variables
# set at their averages):

sd1MCMC <- exp((Xg%*%gammaMCMC + Zg%*%vMCMC)/2)
sd1MCMCOrig <- sd.yObs*sd1MCMC 
credLower <- apply(sd1MCMCOrig,1,quantile,0.025)
credUpper <- apply(sd1MCMCOrig,1,quantile,0.975)
sd1gOrig <- apply(sd1MCMCOrig,1,mean)
plot(x1gOrig,sd1gOrig,type = "n",bty = "l",
     xlab = expression(paste("El Monte temperature (",phantom()^{o},"F)")),
     ylab = "max. 1-hour-average ozone level",
     ylim = ylimValStDev,cex.lab = cex.labVal)
polygon(c(x1gOrig,rev(x1gOrig)),c(credLower,rev(credUpper)),
        col = "palegreen",border = FALSE)
lines(x1gOrig,sd1gOrig,col = "darkgreen",lwd = 2)
rug(x1ObsOrig,col = "dodgerblue")
legend("topleft","standard deviation functions",cex = 1.4)

# Plot second mean function fit (with other variables
# set at their averages):

x1g <- rep(mean(x1Obs),length = ng)
x2g <- seq(min(x2Obs),max(x2Obs),length = ng)
x3g <- rep(mean(x3Obs),length = ng)
Z1g <- outer(x1g,knots1,"-") ; Z1g <- Z1g*(Z1g>0)
Z2g <- outer(x2g,knots2,"-") ; Z2g <- Z2g*(Z2g>0)
Z3g <- outer(x3g,knots3,"-") ; Z3g <- Z3g*(Z3g>0)
Xg <- cbind(rep(1,ng),x1g,x2g,x3g)
Zg <- cbind(Z1g,Z2g,Z3g)
f2MCMC <- Xg%*%betaMCMC + Zg%*%uMCMC
f2MCMCOrig <- sd.yObs*f2MCMC + mean.yObs
credLower <- apply(f2MCMCOrig,1,quantile,0.025)
credUpper <- apply(f2MCMCOrig,1,quantile,0.975)
f2gOrig <- apply(f2MCMCOrig,1,mean)

credLower[credLower<min(ylimValMean)] <- min(ylimValMean)
f2gOrig[f2gOrig<min(ylimValMean)] <- NA

x2gOrig <- sd.x2Obs*x2g + mean.x2Obs
plot(x2gOrig,f2gOrig,type = "n",bty = "l",
     xlab = "Daggett pressure gradient (mm Hg)",
     ylab = "max. 1-hour-average ozone level",ylim = ylimValMean,
     cex.lab = cex.labVal)
polygon(c(x2gOrig,rev(x2gOrig)),c(credLower,rev(credUpper)),
        col = "palegreen",border = FALSE)
lines(x2gOrig,f2gOrig,col = "darkgreen",lwd = 2)
rug(x2ObsOrig,col = "dodgerblue")

# Plot second standard deviation function fit (with other variables
# set at their averages):

sd2MCMC <- exp((Xg%*%gammaMCMC + Zg%*%vMCMC)/2)
sd2MCMCOrig <- sd.yObs*sd2MCMC 
credLower <- apply(sd2MCMCOrig,1,quantile,0.025)
credUpper <- apply(sd2MCMCOrig,1,quantile,0.975)
sd2gOrig <- apply(sd2MCMCOrig,1,mean)
plot(x2gOrig,sd2gOrig,type = "n",bty = "l",
     xlab = "Daggett pressure gradient (mm Hg)",
     ylab = "max. 1-hour-average ozone level",ylim = ylimValStDev,
     cex.lab = cex.labVal)
polygon(c(x2gOrig,rev(x2gOrig)),c(credLower,rev(credUpper)),
        col = "palegreen",border = FALSE)
lines(x2gOrig,sd2gOrig,col = "darkgreen",lwd = 2)
rug(x2ObsOrig,col = "dodgerblue")

# Plot third mean function fit (with other variables
# set at their averages):

x1g <- rep(mean(x1Obs),length = ng)
x2g <- rep(mean(x2Obs),length = ng)
x3g <- seq(min(x3Obs),max(x3Obs),length = ng)
Z1g <- outer(x1g,knots1,"-") ; Z1g <- Z1g*(Z1g>0)
Z2g <- outer(x2g,knots2,"-") ; Z2g <- Z2g*(Z2g>0)
Z3g <- outer(x3g,knots3,"-") ; Z3g <- Z3g*(Z3g>0)
Xg <- cbind(rep(1,ng),x1g,x2g,x3g)
Zg <- cbind(Z1g,Z2g,Z3g)
f3MCMC <- Xg%*%betaMCMC + Zg%*%uMCMC
f3MCMCOrig <- sd.yObs*f3MCMC + mean.yObs
credLower <- apply(f3MCMCOrig,1,quantile,0.025)
credUpper <- apply(f3MCMCOrig,1,quantile,0.975)
f3gOrig <- apply(f3MCMCOrig,1,mean)
x3gOrig <- sd.x3Obs*x3g + mean.x3Obs
plot(x3gOrig,f3gOrig,type = "n",bty = "l",
     xlab = "inversion base height (feet)",
     ylab = "max. 1-hour-average ozone level",ylim = ylimValMean,
     cex.lab = cex.labVal)
polygon(c(x3gOrig,rev(x3gOrig)),c(credLower,rev(credUpper)),
        col = "palegreen",border = FALSE)
lines(x3gOrig,f3gOrig,col = "darkgreen",lwd = 2)
rug(x3ObsOrig,col = "dodgerblue")

# Plot third standard deviation function fit (with other variables
# set at their averages):

sd3MCMC <- exp((Xg%*%gammaMCMC + Zg%*%vMCMC)/2)
sd3MCMCOrig <- sd.yObs*sd3MCMC 
credLower <- apply(sd3MCMCOrig,1,quantile,0.025)
credUpper <- apply(sd3MCMCOrig,1,quantile,0.975)
sd3gOrig <- apply(sd3MCMCOrig,1,mean)
plot(x3gOrig,sd3gOrig,type = "n",bty = "l",
     xlab = "inversion base height (feet)",
     ylab = "max. 1-hour-average ozone level",ylim = ylimValStDev,
     cex.lab = cex.labVal)
polygon(c(x3gOrig,rev(x3gOrig)),c(credLower,rev(credUpper)),
        col = "palegreen",border = FALSE)
lines(x3gOrig,sd3gOrig,col = "darkgreen",lwd = 2)
rug(x3ObsOrig,col = "dodgerblue")

# Do summary plots for unobserved data samples and posterior densities:

yUnobsOrigMCMC <- mean.yObs + sd.yObs*yUnobsMCMC
x1UnobsOrigMCMC <- mean.x1Obs + sd.x1Obs*x1UnobsMCMC
x2UnobsOrigMCMC <- mean.x2Obs + sd.x2Obs*x2UnobsMCMC
x3UnobsOrigMCMC <- mean.x3Obs + sd.x3Obs*x3UnobsMCMC

parms <- list(cbind(yUnobsOrigMCMC[1,],x1UnobsOrigMCMC[1,],
                    x2UnobsOrigMCMC[1,],x3UnobsOrigMCMC[1,]))

parNamesVal <- list(c("1st unobserved","max. 1-hour ave.","ozone level"),
                    c("1st unobserved","El Monte","temperature"),
                    c("1st unobserved","Daggett pressure","gradient"),
                    c("1st unobserved","inversion base","height"))

summMCMC(parms,parNames = parNamesVal,numerSummCex = 1)

########## End of ozoneAna ##########

