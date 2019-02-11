########## R script: proteinMargAddInt.Rs ##########

# For fitting a marginal additive/interaction model to the
# protein data.

# Last changed: 22 AUG 2017

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Set MCMC sample size paramaters:

nWarm <- 1000
nKept <- 1000
nThin <- 1

# Load required packages:

library(HRW)  ;  library(rstan)  

# Extract regression data: 

data(protein)

# Store original versions and their mean and standard deviations:

x1Orig <- protein$BMI
x2Orig <- protein$proteinRecall
yOrig <- protein$proteinBioM
m <- 294 ; n <- 2

mean.x1 <- mean(x1Orig) ; sd.x1 <- sd(x1Orig)
mean.x2 <- mean(x2Orig) ; sd.x2 <- sd(x2Orig)
mean.y <- mean(yOrig) ; sd.y <- sd(yOrig)

# Set regression variables to be standardised versions
# of the protein data frame:

x1 <- (x1Orig - mean.x1)/sd.x1
x2 <- (x2Orig - mean.x2)/sd.x2
x3 <- protein$female
y <- (yOrig - mean.y)/sd.y

# Set up basic variables for the spline component:

a1 <- 1.05*min(x1) - 0.05*max(x1) ; b1 <- 1.05*max(x1) - 0.05*min(x1)
numIntKnots1 <- 23
intKnots1 <- seq(a1,b1,length=(numIntKnots1+2))[-c(1,numIntKnots1+2)]

a2 <- 1.05*min(x2) - 0.05*max(x2) ; b2 <- 1.05*max(x2) - 0.05*min(x2)
numIntKnots2 <- 23
intKnots2 <- seq(a2,b2,length=(numIntKnots2+2))[-c(1,numIntKnots2+2)]

# Obtain the spline component of the Z matrix:

range.x1 <- c(a1,b1) ; range.x2 <- c(a2,b2)
Z1 <- ZOSull(x1,range.x1,intKnots1)
Z2 <- ZOSull(x2,range.x2,intKnots2)

ncZ1 <- ncol(Z1) 
ncZ2 <- ncol(Z2) 

# Form additional Z matrices for interaction terms:

Z10 <- (1-x3)*Z1 ; Z11 <- x3*Z1

# Put regression data in matrix format:

x1Mat <- matrix(NA,m,n) ; x2Mat <- matrix(NA,m,n)
x3Mat <- matrix(NA,m,n) ; yMat <- matrix(NA,m,n)
for (i in 1:m)
{
   x1Mat[i,] <- x1[(i-1)*n+(1:n)]
   x2Mat[i,] <- x2[(i-1)*n+(1:n)]
   x3Mat[i,] <- x3[(i-1)*n+(1:n)]
   yMat[i,] <- y[(i-1)*n+(1:n)]
}

# Set hyperparameters:

sigmaBeta <- 1e5 ; Asigma <- 1e5 ; ASigma <- 1e5

# Specify model in Stan:

margAddIntModModel <- 
'data
{
   int<lower=1> m;             int<lower=1> n;
   int<lower=1> ncZ1;          int<lower=1> ncZ2;  
   real<lower=0> sigmaBeta;    
   real<lower=0> Asigma;       real<lower=0> ASigma;          
   matrix[m,n] x1Mat;          matrix[m,n] x2Mat; 
   matrix[m,n] x3Mat;          matrix[m,n] yMat;
   matrix[m*n,ncZ1] Z1;        matrix[m*n,ncZ2] Z2;
   matrix[m*n,ncZ1] Z10;       matrix[m*n,ncZ1] Z11;
}
parameters
{
   real beta0;              real beta1;
   real beta2;              
   real beta0FvsM;          real beta1FvsM;
   vector[ncZ1] u10;        vector[ncZ1] u11;
   vector[ncZ2] u2; 
   real<lower=0> sigma10;   real<lower=0> sigma11; 
   real<lower=0> sigma2; 
   cov_matrix[n] Sigma;     vector[n] a;
}
transformed parameters
{
   matrix[m,n] meanFunc;
   for (i in 1:m) 
      for (j in 1:n) 
         meanFunc[i,j] =  beta0 + beta1*x1Mat[i,j] + beta2*x2Mat[i,j] 
                                 + beta0FvsM*x3Mat[i,j] 
                                 + beta1FvsM*x1Mat[i,j]*x3Mat[i,j] 
                                 + dot_product(u10,Z10[(i-1)*n+j])
                                 + dot_product(u11,Z11[(i-1)*n+j])
                                 + dot_product(u2,Z2[(i-1)*n+j]);
}
model
{
   vector[n] scaleSigma;

   for (i in 1:m) 
      yMat[i] ~ multi_normal(meanFunc[i],Sigma);	
   	
   u10 ~ normal(0,sigma10);   u11 ~ normal(0,sigma11);
   u2 ~ normal(0,sigma2);

   a ~ inv_gamma(0.5,pow(ASigma,-2));     
   for (j in 1:n) scaleSigma[j] = 4/a[j];
   Sigma ~ inv_wishart(n+1,diag_matrix(scaleSigma));

   beta0 ~ normal(0,sigmaBeta) ;   beta1 ~ normal(0,sigmaBeta);
   beta2 ~ normal(0,sigmaBeta) ;   beta0FvsM ~ normal(0,sigmaBeta);
   beta1FvsM ~ normal(0,sigmaBeta) ;
   sigma10 ~ cauchy(0,Asigma);     sigma11 ~ cauchy(0,Asigma);
   sigma2 ~ cauchy(0,Asigma);      

}'

# Set up input data for Stan:

allData <- list(m = m,n = n,ncZ1 = ncZ1,ncZ2 = ncZ2,x1Mat = x1Mat,x2Mat = x2Mat,
                x3Mat = x3Mat,yMat = yMat,Z1 = Z1,Z10 = Z10,Z11 = Z11,Z2 = Z2, 
                sigmaBeta = sigmaBeta,Asigma = Asigma,ASigma = ASigma)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = margAddIntModModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = margAddIntModModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 100,
                 fit = stanCompilObj)

# Save and extract relevant MCMC samples:

beta0MCMC <- as.vector(extract(stanObj,"beta0",permuted = FALSE))
beta1MCMC <- as.vector(extract(stanObj,"beta1",permuted = FALSE))
beta2MCMC <- as.vector(extract(stanObj,"beta2",permuted = FALSE))
beta0FvsMMCMC <- as.vector(extract(stanObj,"beta0FvsM",permuted = FALSE))
beta1FvsMMCMC <- as.vector(extract(stanObj,"beta1FvsM",permuted = FALSE))
betaMCMC <- rbind(beta0MCMC,beta1MCMC,beta0FvsMMCMC,beta1FvsMMCMC,beta2MCMC)

u10MCMC <- NULL    ;  u11MCMC <- NULL ; u2MCMC <- NULL  

for (k in 1:ncZ1)
{
   charVar <- paste("u10[",as.character(k),"]",sep = "") 
   u10MCMC <- rbind(u10MCMC,extract(stanObj,charVar,permuted = FALSE))
   charVar <- paste("u11[",as.character(k),"]",sep="") 
   u11MCMC <- rbind(u11MCMC,extract(stanObj,charVar,permuted = FALSE))
}
for (k in 1:ncZ2)
{
   charVar <- paste("u2[",as.character(k),"]",sep="") 
   u2MCMC <- rbind(u2MCMC,extract(stanObj,charVar,permuted = FALSE))
}

sigma10MCMC <- as.vector(extract(stanObj,"sigma10",permuted = FALSE))
sigma11MCMC <- as.vector(extract(stanObj,"sigma11",permuted = FALSE))
sigma2MCMC <- as.vector(extract(stanObj,"sigma2",permuted = FALSE))

Sigma11MCMC <- extract(stanObj,"Sigma[1,1]",permuted = FALSE)
Sigma12MCMC <- extract(stanObj,"Sigma[1,2]",permuted = FALSE)
Sigma22MCMC <- extract(stanObj,"Sigma[2,2]",permuted = FALSE)

# Transform MCMC samples to correspond to original units:

sigma10MCMCOrig <- sigma10MCMC*sd.y
sigma11MCMCOrig <- sigma11MCMC*sd.y
sigma2MCMCOrig <- sigma2MCMC*sd.y

Sigma11MCMCOrig <-  Sigma11MCMC*(sd.y^2)
Sigma22MCMCOrig <-  Sigma22MCMC*(sd.y^2)
Sigma12MCMCOrig <-  Sigma12MCMC*(sd.y^2)

# Plot curve estimates:

par(mfrow=c(1,3),mai=c(1.1,1.2,0.1,0.2),mgp=c(6,2,0))
cex.labVal <- 3.2  ; cex.axisVal <- 2.8
ng <- 101 ; ylimVal <- c(5.2,6.6)

shade <- FALSE ; colourVersion <- FALSE

meanCol <- "darkgreen" ; seCol <- "palegreen"
rugCol <- "dodgerblue"
      
# First function plots:

x1g <- seq(range.x1[1],range.x1[2],length=ng)
x2g <- rep(mean(x2),length=ng)
Z1g <- ZOSull(x1g,range.x1,intKnots1)
Z2g <- ZOSull(x2g,range.x2,intKnots2)

# Plot for female=0:

X0g <- cbind(rep(1,ng),x1g,rep(0,ng),rep(0,ng),x2g)
f10MCMC <- X0g%*%betaMCMC + Z1g%*%u10MCMC + Z2g%*%u2MCMC
f10MCMCorig <- mean.y + sd.y*f10MCMC
f10gOrig <- apply(f10MCMCorig,1,mean)
credLower10Orig <- apply(f10MCMCorig,1,quantile,0.025)
credUpper10Orig <- apply(f10MCMCorig,1,quantile,0.975)
x1gOrig <- mean.x1 + sd.x1*x1g

plot(x1gOrig,f10gOrig,type="n",ylim=ylimVal,
     xlab="body mass index (males)",cex.axis=cex.axisVal,
     ylab="mean log(protein biomarker)",cex.lab=cex.labVal,
     bty="l")

polygon(c(x1gOrig,rev(x1gOrig)),c(credLower10Orig,rev(credUpper10Orig)),
        col=seCol,border=FALSE)

lines(x1gOrig,f10gOrig,lwd=2,col=meanCol)
rug(x1Orig,quiet=TRUE,col=rugCol)

# Plot for female=1:

X1g <- cbind(rep(1,ng),x1g,rep(1,ng),rep(1,ng),x2g)
f11MCMC <- X1g%*%betaMCMC + Z1g%*%u11MCMC + Z2g%*%u2MCMC
f11MCMCorig <-  mean.y + sd.y*f11MCMC
f11gOrig <- apply(f11MCMCorig,1,mean)
credLower11Orig <- apply(f11MCMCorig,1,quantile,0.025)
credUpper11Orig <- apply(f11MCMCorig,1,quantile,0.975)
x1gOrig <- mean.x1 + sd.x1*x1g

plot(x1gOrig,f11gOrig,type="n",ylim=ylimVal,cex.axis=cex.axisVal,
     cex.lab=cex.labVal,
     xlab="body mass index (females)",ylab="mean log(protein biomarker)",
     bty="l")

polygon(c(x1gOrig,rev(x1gOrig)),c(credLower11Orig,rev(credUpper11Orig)),
                 col=seCol,border=FALSE)

lines(x1gOrig,f11gOrig,lwd=2,col=meanCol)
rug(x1Orig,quiet=TRUE,col=rugCol)

# Plot second function fit (with other variable
# set at its average):

x1g <- rep(mean(x1),length=ng)
x2g <- seq(range.x2[1],range.x2[2],length=ng)
Xg <- cbind(rep(1,ng),x1g,rep(1,ng),x1g,x2g)
Z1g <- ZOSull(x1g,range.x1,intKnots1)
Z2g <- ZOSull(x2g,range.x2,intKnots2)

f2MCMC <- Xg%*%betaMCMC + Z1g%*%u11MCMC + Z2g%*%u2MCMC
f2MCMCorig <-  mean.y + sd.y*f2MCMC
f2gOrig <- apply(f2MCMCorig,1,mean)
credLower2Orig <- apply(f2MCMCorig,1,quantile,0.025)
credUpper2Orig <- apply(f2MCMCorig,1,quantile,0.975)
x2gOrig <- mean.x2 + sd.x2*x2g

plot(x2gOrig,f2gOrig,type="n",ylim=ylimVal,cex.axis=cex.axisVal,
     cex.lab=cex.labVal,xlab="protein recall",ylab="mean log(protein biomarker)",
     bty="l")

polygon(c(x2gOrig,rev(x2gOrig)),c(credLower2Orig,rev(credUpper2Orig)),
        col=seCol,border=FALSE)

lines(x2gOrig,credLower2Orig,lwd=2,lty=2,col=seCol)
lines(x2gOrig,credUpper2Orig,lwd=2,lty=2,col=seCol)

lines(x2gOrig,f2gOrig,lwd=2,col=meanCol)
rug(x2Orig,quiet=TRUE,col=rugCol)

# Do MCMC summary plot:

indQ2 <- length(x1gOrig[x1gOrig<quantile(x1Orig,0.50)])
fhat10OrigQ2 <- f10MCMCorig[indQ2,]
fhat11OrigQ2 <- f11MCMCorig[indQ2,]

indQ2 <- length(x2gOrig[x2gOrig<quantile(x2Orig,0.50)])
fhat2OrigQ2 <- f2MCMCorig[indQ2,]

parmsMCMC <- list(cbind(fhat10OrigQ2,fhat11OrigQ2,fhat2OrigQ2,
                        Sigma11MCMCOrig,Sigma12MCMCOrig,Sigma22MCMCOrig))
parNamesVec <- list(c("mean function","at median","of BMI (males)"),
                    c("mean function","at median","of BMI (females)"),
                    c("mean function","at median","of protein recall"),
                    expression(Sigma[11]),expression(Sigma[12]),
                    expression(Sigma[22]))

summMCMC(parmsMCMC,parNames=parNamesVec,KDEvertLine=TRUE)

########## End of proteinMargAddInt ##########


