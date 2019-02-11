########## R script: margAddMod ##########

# For fitting a marginal additive model to simulated data.

# Last changed: 22 AUG 2017

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Set MCMC sample size paramaters:

nWarm <- 1000
nKept <- 1000
nThin <- 1

# Load required packages:

library(HRW)   ;   library(rstan)   ;   library(lattice)

# Define a matrix square-root function:

matrixSqrt <- function(A)
{
   sva <- svd(A)
   if (min(sva$d) >= 0)
      Asqrt <- t(sva$v%*%(t(sva$u)*sqrt(sva$d)))
   else
      stop("Matrix square root is not defined")
   return(Asqrt)
}

# Simulate some data:

set.seed(44601) ; m <- 200 ; n <- 5

f1 <- function(x)
    return(sin(2*pi*(x^2-0.1)))

f2 <- function(x)
    return(sin(3*pi*(0.05-x)))

corrMat <- toeplitz(c(0.97,0.83,0.77,0.41,0.18))
SigmaTrue <- corrMat
sqrtSigmaTrue <- matrixSqrt(SigmaTrue)

x1 <- NULL ; x2 <- NULL 
idNum <- NULL

x1stts <- runif(m,0,1/n) 
x2stts <- runif(m,0,1/n) 

epsVals <- NULL
for (i in 1:m)
{
   idNum <- c(idNum,rep(i,n))
   x1 <- c(x1,seq(x1stts[i],by = 1/n,length = n))
   x2 <- c(x2,seq(x2stts[i],by = 1/n,length = n))
   epsVals <- c(epsVals,as.vector(sqrtSigmaTrue%*%rnorm(n)))
}

y <- f1(x1) + f2(x2) + epsVals

# Convert y to matrix format:

x1Mat <- matrix(NA,m,n) ; x2Mat <- matrix(NA,m,n)
yMat <- matrix(NA,m,n)
for (i in 1:m)
{
   x1Mat[i,] <- x1[(i-1)*n+(1:n)]
   x2Mat[i,] <- x2[(i-1)*n+(1:n)]
   yMat[i,] <- y[(i-1)*n+(1:n)]
}

# Set up basic variables for the spline components:

numIntKnots1 <- 15 ; range.x1 <- c(0,1)
intKnots1 <- seq(range.x1[1],range.x1[2],
                 length = (numIntKnots1+2))[-c(1,numIntKnots1+2)]

numIntKnots2 <- 20 ; range.x2 <- c(0,1)
intKnots2 <- seq(range.x2[1],range.x2[2],
                 length = (numIntKnots2+2))[-c(1,numIntKnots2+2)]

# Obtain the spline component of the Z matrix:

Z1 <- ZOSull(x1,range.x1,intKnots1)
Z2 <- ZOSull(x2,range.x2,intKnots2)

ncZ1 <- ncol(Z1) 
ncZ2 <- ncol(Z2) 

# Set hyperparameters:

sigmaBeta <- 1e10 ; Asigma <- 1e10 ; ASigma <- 1e10

# Specify model in Stan:

margAddModModel <- 
'data
{
   int<lower = 1> m;             int<lower=1> n;
   int<lower=1> ncZ1;          int<lower=1> ncZ2;  
   real<lower=0> sigmaBeta;    
   real<lower=0> Asigma;       real<lower=0> ASigma;          
   matrix[m,n] x1Mat;          matrix[m,n] x2Mat; 
   matrix[m,n] yMat;
   matrix[m*n,ncZ1] Z1;        matrix[m*n,ncZ2] Z2;
}
parameters
{
   real beta0;              real beta1;
   real beta2;
   vector[ncZ1] u1;         vector[ncZ2] u2; 
   real<lower=0> sigma1;    real<lower=0> sigma2; 
   cov_matrix[n] Sigma;     vector[n] a;
}
transformed parameters
{
   matrix[m,n] meanFunc;
   for (i in 1:m) 
      for (j in 1:n) 
         meanFunc[i,j] = beta0 + beta1*x1Mat[i,j] + dot_product(u1,Z1[(i-1)*n+j])
                               + beta2*x2Mat[i,j] + dot_product(u2,Z2[(i-1)*n+j]);
}
model
{
   matrix[n,n] rateSigma;

   for (i in 1:m) 
      yMat[i] ~ multi_normal(meanFunc[i],Sigma);	
   	
   u1 ~ normal(0,sigma1);
   u2 ~ normal(0,sigma2);

   rateSigma = rep_matrix(0,n,n);
   for (j in 1:n)
   {
      a[j] ~ inv_gamma(0.5,pow(ASigma,-2));     
      rateSigma[j,j] = 4/a[j];
   }

   Sigma ~ inv_wishart(n+1,rateSigma);

   beta0  ~ normal(0,sigmaBeta) ; beta1 ~ normal(0,sigmaBeta);
   beta2 ~ normal(0,sigmaBeta);
   sigma1 ~ cauchy(0,Asigma);     sigma2 ~ cauchy(0,Asigma);
}'

# Obtain Markov Chain Monte Carlo (MCMC) samples:

allData <- list(m = m,n = n,ncZ1 = ncZ1,ncZ2 = ncZ2,x1Mat = x1Mat,x2Mat = x2Mat,
                yMat = yMat,Z1 = Z1,Z2 = Z2,sigmaBeta = sigmaBeta,
                Asigma = Asigma,ASigma = ASigma)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = margAddModModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = margAddModModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 100,
                 fit = stanCompilObj)

# Save and extract relevant MCMC samples:

beta0MCMC <- extract(stanObj,"beta0",permuted = FALSE)
beta1MCMC <- extract(stanObj,"beta1",permuted = FALSE)
beta2MCMC <- extract(stanObj,"beta2",permuted = FALSE)

betaMCMC <- rbind(beta0MCMC,beta1MCMC,beta2MCMC)

u1MCMC <- NULL
for (k in 1:ncZ1)
{
   charVar <- paste("u1[",as.character(k),"]",sep = "") 
   u1MCMC <- rbind(u1MCMC,extract(stanObj,charVar,permuted = FALSE))
}

u2MCMC <- NULL
for (k in 1:ncZ2)
{
   charVar <- paste("u2[",as.character(k),"]",sep = "") 
   u2MCMC <- rbind(u2MCMC,extract(stanObj,charVar,permuted = FALSE))
}

sigma1MCMC <- extract(stanObj,"sigma1",permuted = FALSE)
sigma2MCMC <- extract(stanObj,"sigma2",permuted = FALSE)
Sigma11MCMC <- extract(stanObj,"Sigma[1,1]",permuted = FALSE)
Sigma22MCMC <- extract(stanObj,"Sigma[2,2]",permuted = FALSE)
Sigma33MCMC <- extract(stanObj,"Sigma[3,3]",permuted = FALSE)
Sigma44MCMC <- extract(stanObj,"Sigma[4,4]",permuted = FALSE)
Sigma55MCMC <- extract(stanObj,"Sigma[5,5]",permuted = FALSE)
Sigma12MCMC <- extract(stanObj,"Sigma[1,2]",permuted = FALSE)
Sigma13MCMC <- extract(stanObj,"Sigma[1,3]",permuted = FALSE)
Sigma14MCMC <- extract(stanObj,"Sigma[1,4]",permuted = FALSE)
Sigma15MCMC <- extract(stanObj,"Sigma[1,5]",permuted = FALSE)
Sigma23MCMC <- extract(stanObj,"Sigma[2,3]",permuted = FALSE)
Sigma24MCMC <- extract(stanObj,"Sigma[2,4]",permuted = FALSE)
Sigma25MCMC <- extract(stanObj,"Sigma[2,5]",permuted = FALSE)
Sigma34MCMC <- extract(stanObj,"Sigma[3,4]",permuted = FALSE)
Sigma35MCMC <- extract(stanObj,"Sigma[3,5]",permuted = FALSE)
Sigma45MCMC <- extract(stanObj,"Sigma[4,5]",permuted = FALSE)

# Plot the MCMC results, starting with diagonal entries
# of the Sigma matrix:

parmsMCMC <- list(cbind(Sigma11MCMC,Sigma22MCMC,Sigma33MCMC,Sigma44MCMC,Sigma55MCMC))
parNamesVec <- list(expression(Sigma[11]),expression(Sigma[22]),
                    expression(Sigma[33]),expression(Sigma[44]),expression(Sigma[55]))

summMCMC(parmsMCMC,parNames=parNamesVec,KDEvertLine=FALSE,
         addTruthToKDE=c(SigmaTrue[1,1],SigmaTrue[2,2],SigmaTrue[3,3],   
                         SigmaTrue[4,4],SigmaTrue[5,5]))

# Off-diagonal entries:

parmsMCMC <- list(cbind(Sigma11MCMC,Sigma12MCMC,Sigma13MCMC,Sigma14MCMC,
                        Sigma23MCMC,Sigma24MCMC,Sigma25MCMC,Sigma34MCMC,
                        Sigma35MCMC,Sigma45MCMC))

parNamesVec <- list(expression(Sigma[12]),expression(Sigma[13]),
                    expression(Sigma[14]),expression(Sigma[15]),
                    expression(Sigma[23]),expression(Sigma[24]),
                    expression(Sigma[25]),expression(Sigma[34]),
                    expression(Sigma[35]),expression(Sigma[45]))

summMCMC(parmsMCMC,parNames=parNamesVec,KDEvertLine = FALSE,
         addTruthToKDE = c(SigmaTrue[1,2],SigmaTrue[1,3],SigmaTrue[1,4],   
                         SigmaTrue[1,5],SigmaTrue[2,3],SigmaTrue[2,4],
                         SigmaTrue[2,5],SigmaTrue[3,4],SigmaTrue[3,5],
                         SigmaTrue[4,5]))

# Plot first function fit (with other variable
# set at its average):

par(mfrow = c(1,2),mai = c(0.9,0.76,0.1,0.1))

ylimVal <- c(-0.8,2.2)  ;   ng <- 201
meanCol <- "darkgreen"  ;  seCol <- "palegreen"  ;  truthCol <- "indianred3"
dataCol <- "dodgerblue"

# First function:

x1g <- seq(range.x1[1],range.x1[2],length = ng)
x2g <- rep(mean(x2),length = ng)

Z1g <- ZOSull(x1g,range.x1,intKnots1)
Z2g <- ZOSull(x2g,range.x2,intKnots2)

Xg <- cbind(rep(1,ng),x1g,x2g)
Zg <- cbind(Z1g,Z2g)
f1MCMC <- Xg%*%betaMCMC + Z1g%*%u1MCMC + Z2g%*%u2MCMC
credLower <- apply(f1MCMC,1,quantile,0.025)
credUpper <- apply(f1MCMC,1,quantile,0.975)
f1g <- apply(f1MCMC,1,mean)
plot(x1g,f1g,type = "n",xlab = expression(x[1]),ylab = "mean y",ylim = ylimVal,bty = "l")
polygon(c(x1g,rev(x1g)),c(credLower,rev(credUpper)),col = seCol,border = FALSE)
rug(x1,col = "dodgerblue")
lines(x1g,f1g,lwd = 2,col = meanCol)
f1gTrue <- f1(x1g) + f2(x2g)
lines(x1g,f1gTrue,col = truthCol)

legend("topleft",legend = c(expression(paste("true ",f[1])),
       expression(paste("estimated ",f[1]))),
       lwd = c(1,2),col = c("red","forestgreen"))

# Plot second function fit (with other variable set at its average):

x1g <- rep(mean(x1),length = ng)
x2g <- seq(range.x2[1],range.x2[2],length = ng)

Z1g <- ZOSull(x1g,range.x1,intKnots1)
Z2g <- ZOSull(x2g,range.x2,intKnots2)

Xg <- cbind(rep(1,ng),x1g,x2g)
Zg <- cbind(Z1g,Z2g)
f2MCMC <- Xg%*%betaMCMC + Z1g%*%u1MCMC + Z2g%*%u2MCMC
credLower <- apply(f2MCMC,1,quantile,0.025)
credUpper <- apply(f2MCMC,1,quantile,0.975)
f2g <- apply(f2MCMC,1,mean)
plot(x2g,f2g,type = "n",xlab = expression(x[2]),ylab = "mean y",ylim = ylimVal,bty = "l")
polygon(c(x2g,rev(x2g)),c(credLower,rev(credUpper)),col = seCol,border = FALSE)
rug(x2,col = dataCol)
lines(x2g,f2g,lwd = 2,col = meanCol)
f2gTrue <- f1(x1g) + f2(x2g)
lines(x2g,f2gTrue,col = truthCol)

############ End of margAddMod ############
