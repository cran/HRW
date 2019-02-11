########## R script: maleGrowthIndianaBayes ##########

# For performing a Bayesian group-specific curve analysis
# of the Indiana growth data for male subjects.

# Last changed: 29 AUG 2018

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Set MCMC sample size paramaters:

nWarm <- 1000
nKept <- 1000
nThin <- 1
   
# Load required packages:

library(HRW)   ;   library(rstan)   ;  library(lattice) 

# Load in the data:

data(growthIndiana)

# Extract data on males:

growthIndianaMales <- growthIndiana[growthIndiana$male == 1,]

# Create relevant variables: 

height <- growthIndianaMales$height
age <- growthIndianaMales$age
idnumOrig <- growthIndianaMales$idnum
typeIsB <- growthIndianaMales$black

# Create new (ordered) ID numbers:

idnum <- rep(NA,length(idnumOrig))
uqID <- unique(idnumOrig)
for (i in 1:length(uqID))
      idnum[idnumOrig == uqID[i]] <- i
  
# Use xOrig and yOrig to denote the original continuous
# predictor and response variables:

xOrig <- age ;  yOrig <- height

# Save mean and standard deviation information:

mean.x <- mean(xOrig)   ;   sd.x <- sd(xOrig)
mean.y <- mean(yOrig)   ;   sd.y <- sd(yOrig)

# Store the standardised versions of the predictor
# and response variables in x and y:

x <- (xOrig - mean.x)/sd.x
y <- (yOrig - mean.y)/sd.y

# Set up X and Z matrix for mean curve:

numObs <- length(x)
numGrp <- length(unique(idnum))

Xbase <- cbind(rep(1,numObs),x)
XB <- typeIsB*Xbase

numIntKnots <- 20
intKnots <- quantile(unique(x),
                     seq(0,1,length=numIntKnots+2))[-c(1,numIntKnots+2)]

range.x <- c(1.01*min(x) - 0.01*max(x),1.01*max(x)-0.01*min(x))
 
Z <- ZOSull(x,intKnots=intKnots,range.x=range.x)
ZA <- (1-typeIsB)*Z
ZB <- typeIsB*Z

# Set up Z matrix for group specific curves:

numIntKnotsGrp <- 10
intKnotsGrp <- quantile(unique(x),
        seq(0,1,length=numIntKnotsGrp+2))[-c(1,numIntKnotsGrp+2)]

Zgrp <- ZOSull(x,intKnots=intKnotsGrp,range.x=range.x)

# Specify hyperparameters:

sigmaBeta <- 1e5  ;  AU <- 1e5
AuGbl <- 1e5      ;  AuGrp <- 1e5
Aeps <- 1e5 

# Let up dimension variables:
   
ncZ <- ncol(Z) ; ncZgrp <- ncol(Zgrp)

# Set up model in Stan:

normGrpSpecContrModel <-
'data
{
   int<lower=1> numObs;             int<lower=1> numGrp;
   int<lower=1> ncZ;                int<lower=1> ncZgrp;
   real<lower=0> sigmaBeta;         real<lower=0> AU;
   real<lower=0> AuGbl;             real<lower=0> AuGrp;
   real<lower=0> Aeps;              matrix[numObs,2] Xbase;          
   matrix[numObs,2] XB;             real y[numObs];                  
   int<lower=1> idnum[numObs];      matrix[numObs,ncZ] ZA;        
   matrix[numObs,ncZ] ZB;           matrix[numObs,ncZgrp] Zgrp;
}
transformed data
{
   vector[2] zeroVec;           zeroVec = rep_vector(0,2);
}
parameters
{
   vector[2] beta;               vector[2] betaBvsA;
   vector[2] U[numGrp];          vector[ncZ] uGblA;
   vector[ncZ] uGblB;            matrix[numGrp,ncZgrp] uGrp;
   cov_matrix[2] SigmaGrp;       real<lower=0> siguGblA;
   real<lower=0> siguGblB;       real<lower=0> siguGrp;
   real<lower=0> sigEps;         vector[2] a;
}
transformed parameters
{
   vector[numObs] meanVec;
   meanVec = Xbase*beta + XB*betaBvsA + ZA*uGblA + ZB*uGblB   
             + to_vector(U[idnum,1]) + to_vector(U[idnum,2]).*col(Xbase,2)
             + rows_dot_product(uGrp[idnum],Zgrp);
}
model
{
   vector[2] scaleSigmaGrp;

   y ~ normal(meanVec,sigEps);

   U ~ multi_normal(zeroVec,SigmaGrp);
   uGblA ~ normal(0,siguGblA);  uGblB ~ normal(0,siguGblB);     
   to_vector(uGrp) ~ normal(0,siguGrp);
      
   a ~ inv_gamma(0.5,pow(AU,-2)); 
   for (k in 1:2) scaleSigmaGrp[k] = 4/a[k];
   SigmaGrp ~ inv_wishart(3,diag_matrix(scaleSigmaGrp));

   beta ~ normal(0,sigmaBeta);   betaBvsA ~ normal(0,sigmaBeta);
   siguGblA ~ cauchy(0,AuGbl);   siguGblB ~ cauchy(0,AuGbl);
   siguGrp ~ cauchy(0,AuGrp);    sigEps ~ cauchy(0,Aeps);
}'

	  
# Set up data list:

allData <- list(numObs=numObs,numGrp=numGrp,ncZ=ncZ,
                ncZgrp=ncZgrp,AU=AU,AuGbl=AuGbl,
                AuGrp=AuGrp,Aeps=Aeps,sigmaBeta=sigmaBeta,
                y=y,idnum=idnum,ZA=ZA,ZB=ZB,
     	        Zgrp=Zgrp,Xbase=Xbase,XB=XB)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code=normGrpSpecContrModel,data=allData,
                         iter=1,chains=1)

# Perform MCMC:

stanObj <-  stan(model_code=normGrpSpecContrModel,data=allData,warmup=nWarm,
                 iter=(nWarm + nKept),chains=1,thin=nThin,refresh=10,
                 fit=stanCompilObj)

# Save and extract relevant MCMC samples:

betaMCMC <- NULL  ;  betaBvsAMCMC <- NULL   
for (i in 1:2)
{
   charVar1 <- paste("beta[",as.character(i),"]",sep="") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar1,permuted=FALSE))
   charVar2 <- paste("betaBvsA[",as.character(i),"]",sep="") 
   betaBvsAMCMC <- rbind(betaBvsAMCMC,extract(stanObj,charVar2,permuted=FALSE))
}

uGblAMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("uGblA[",as.character(k),"]",sep="") 
   uGblAMCMC <- rbind(uGblAMCMC,extract(stanObj,charVar,permuted=FALSE))
}
uGblBMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("uGblB[",as.character(k),"]",sep="") 
   uGblBMCMC <- rbind(uGblBMCMC,extract(stanObj,charVar,permuted=FALSE))
}
siguGblAMCMC <- as.vector(extract(stanObj,"siguGblA",permuted=FALSE))
siguGblBMCMC <- as.vector(extract(stanObj,"siguGblB",permuted=FALSE))
siguGrpMCMC <- as.vector(extract(stanObj,"siguGrp",permuted=FALSE))
sigEpsMCMC <- as.vector(extract(stanObj,"sigEps",permuted=FALSE))

UMCMC <- NULL
for (i in 1:numGrp)
   for (j in 1:2)
   {
      charVar <- paste("U[",as.character(i),",",as.character(j),"]",sep="")
      UMCMC <- rbind(UMCMC,as.vector(extract(stanObj,charVar,permuted=FALSE)))
   }

uGrpMCMC <- NULL
for (i in 1:numGrp)
   for (j in 1:ncZgrp)
   {
      charVar <- paste("uGrp[",as.character(i),",",as.character(j),"]",sep="")
      uGrpMCMC <- rbind(uGrpMCMC,as.vector(extract(stanObj,charVar,permuted=FALSE)))
   }

# Do lattice-type graphic showing fitted group-specific curves:

AptCol <- "dodgerblue"   ;   AlnCol <- "blue"
BptCol <- "deeppink"     ;   BlnCol <- "maroon"
lwdVal <- 1.5

ng <- 101
xg <- seq(min(x),max(x),length=ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,intKnots=intKnots,range.x=range.x)
Zggrp <- ZOSull(xg,intKnots=intKnotsGrp,range.x=range.x)

fAgMCMC <- Xg%*%betaMCMC + Zg%*%uGblAMCMC
fBgMCMC <- Xg%*%(betaMCMC + betaBvsAMCMC) + Zg%*%uGblBMCMC

meanCrvs <- vector("list",numGrp)
credLows <- vector("list",numGrp)
credUpps <- vector("list",numGrp)
ggMCMC <- vector("list",numGrp)
for (i in 1:numGrp)
{
   # Determine MCMC sample of group-specific deviation
   # grids for group i:

   ggMCMC[[i]] <- (Xg%*%rbind(UMCMC[(2*i-1),],UMCMC[(2*i),])  
                  + Zggrp%*%uGrpMCMC[(i-1)*ncZgrp+1:ncZgrp,])

   # Determine if the current group is Type A or Type B:
  
   indsi <- (1:numObs)[idnum == i]
   currIsB <- typeIsB[indsi[1]]
  
   if (!currIsB)
      currMCMC <-  fAgMCMC + ggMCMC[[i]]

   if (currIsB)
      currMCMC <-  fBgMCMC + ggMCMC[[i]]
   
   meanCrvs[[i]] <- apply(currMCMC,1,mean)
   credLows[[i]] <- apply(currMCMC,1,quantile,0.025)
   credUpps[[i]] <- apply(currMCMC,1,quantile,0.975)
}

# Convert sigEpsMCMC to the original units:

sigEpsMCMCorig <-  sd.y*sigEpsMCMC

# Convert fAgMCMC, fBgMCMC to the original units:

fAgMCMCorig <-  mean.y + sd.y*fAgMCMC
fBgMCMCorig <-  mean.y + sd.y*fBgMCMC

# Convert the group-specific curves to the orginal units:

meanCrvsOrig <- vector("list",numGrp)
credLowsOrig <- vector("list",numGrp)
credUppsOrig <- vector("list",numGrp)

for (i in 1:numGrp)
{
   meanCrvsOrig[[i]] <- mean.y + sd.y*meanCrvs[[i]]
   credLowsOrig[[i]] <- mean.y + sd.y*credLows[[i]]
   credUppsOrig[[i]] <- mean.y + sd.y*credUpps[[i]]
}  
  
# Convert the horizontal plotting grid to the original units:

xgOrig <- mean.x + sd.x*xg

allDataDF <- data.frame(xOrig,yOrig,idnum,typeIsB)

cex.labVal <- 1.35
colourVec <- c("dodgerblue","deeppink","blue","darkmagenta","lightblue","pink")

figFit <- xyplot(yOrig~xOrig|idnum,groups = idnum,data = allDataDF,
                 layout = c(9,13),xlab = list(label = "age (years)",cex = cex.labVal),
                 ylab = list(label = "height (centimetres)",cex = cex.labVal),
                 cex.lab = 4,strip = FALSE,as.table = TRUE,
                 key = list(title = "",
                      columns = 2, 
                      lines = list(lty = rep(1,2),col = colourVec[3:4]),
                      text = list(c("black","white"),cex = 1.55)),
                 panel = function(x,y,subscripts,groups)
                 {
                    panel.grid()
                    iGrp <- panel.number()
                    colourInd <- 2 - typeIsB[subscripts[1]] 
                    panel.superpose(x,y,subscripts,groups,
                                    type = "b",col = colourVec[colourInd],
                                    pch = 1,cex = 0.5)                   
                    panel.polygon(c(xgOrig,rev(xgOrig)),
                                  c(credLowsOrig[[iGrp]],
                                  rev(credUppsOrig[[iGrp]])),
                                  col = colourVec[colourInd+4],border = FALSE)
                    panel.xyplot(xgOrig,meanCrvsOrig[[iGrp]],
                                 col = colourVec[colourInd+2],
                                 type = "l",lwd = lwdVal)
                    
                 })
print(figFit)

# Do contrast function plot:

ContgMCMCorig <- fBgMCMCorig - fAgMCMCorig
ContgMeanOrig <- apply(ContgMCMCorig,1,mean)
ContgLowerOrig <- apply(ContgMCMCorig,1,quantile,0.025)
ContgUpperOrig <- apply(ContgMCMCorig,1,quantile,0.975)

cex.labVal <- 2 ; cex.axisVal <- 1.5
ylimVal <- range(c(ContgLowerOrig,ContgUpperOrig))
par(mai = c(1,0.95,0.1,0.1))
plot(xgOrig,ContgMeanOrig,bty = "l",type = "n",xlab = "age (years)",
     ylab = "mean difference in height (centimetres)",
     ylim = ylimVal,cex.lab = cex.labVal,cex.axis = cex.axisVal)

polygon(c(xgOrig,rev(xgOrig)),c(ContgLowerOrig,rev(ContgUpperOrig)),
         col = "palegreen",border = FALSE)
lines(xgOrig,ContgMeanOrig,col = "darkgreen",lwd = 2)
abline(0,0,col = "slateblue")
nSub <- 1000
subInds <- sample(1:length(x),nSub,replace = FALSE)
rug(xOrig[subInds],col = "dodgerblue")

# Do MCMC summary plot to check convergence:

indQ1 <- length(xgOrig[xgOrig<quantile(xOrig,0.25)])
indQ2 <- length(xgOrig[xgOrig<quantile(xOrig,0.50)])
indQ3 <- length(xgOrig[xgOrig<quantile(xOrig,0.75)])
ContQ1MCMCorig <- ContgMCMCorig[indQ1,]
ContQ2MCMCorig <- ContgMCMCorig[indQ2,]
ContQ3MCMCorig <- ContgMCMCorig[indQ3,]

MCMClist <- list(cbind(sigEpsMCMCorig,ContQ1MCMCorig,ContQ2MCMCorig,ContQ3MCMCorig))
parNamesVal <- list(expression(sigma[epsilon]),
                    c("contrast function","at 1st quantile","of age"),
                    c("contrast function","at 2nd quantile","of age"),
                    c("contrast function","at 3rd quantile","of age"))
summMCMC(MCMClist,parNames = parNamesVal)

############ End of maleGrowthIndianaBayes ###########








