########## R script: coralAna ##########

# For performing a group-specific curve analysis
# of the French Polynesia coral data.

# Last changed: 07 JUN 2018

# Set flag for code compilation (needed if 
# running script first time in current session) :

compileCode <- TRUE

# Set MCMC sample size paramaters:

nWarm <- 1000
nKept <- 1000
nThin <- 1

# Load required packages and functions:

library(lattice)   ;   library(rstan)
library(HRW)

# Set flag for truncation of the x-axis data at 4:

truncatexAt4 <- TRUE

# Load in the data:

data(coral)

# Restrict to Pocillopora and Porites coral:

Ainds <- (1:nrow(coral))[coral$taxon == "POC"]
Binds <- (1:nrow(coral))[coral$taxon == "POR"]
PocPorData <- coral[union(Ainds,Binds),]
PocPorData <- PocPorData[order(PocPorData$siteDepthPeriod),]

# Truncate high (sparse) x values if requested:

if (truncatexAt4)
{
   indsOmit <- (1:nrow(PocPorData))[PocPorData$logInitialSizePlusOne>=4]
   PocPorData <- PocPorData[-indsOmit,]
}

# Create relevant variables: 

taxon <- as.character(PocPorData$taxon)
logInitialSizePlusOne <- PocPorData$logInitialSizePlusOne
died <- PocPorData$died
siteDepthPeriod <- as.character(PocPorData$siteDepthPeriod)

# Create new (ordered) ID numbers:

idnum <- rep(NA,length(siteDepthPeriod))
uqID <- unique(siteDepthPeriod)
for (i in 1:length(uqID))
      idnum[siteDepthPeriod == uqID[i]] <- i
  
# Use x and y to denote the continuous predictor
# and response. Note that, for these data, the 
# predictor values are between 0 and 4 so it is 
# reasonable to work with the original units.

x <- logInitialSizePlusOne
y <- 1 - died

# Obtain taxon indicator variable:

typeIsB <- as.numeric(taxon == "POR") 

# Store dimension variables:

numObs <- length(x)
numGrp <- length(unique(idnum))

# Set up X and Z matrix for mean curve:

Xbase <- cbind(rep(1,numObs),x)
XB <- typeIsB*Xbase
XA <- (1-typeIsB*Xbase)

numIntKnots <- 20
intKnots <- quantile(unique(x),
                     seq(0,1,length = numIntKnots+2))[-c(1,numIntKnots+2)]

range.x <- c(1.01*min(x) - 0.01*max(x),1.01*max(x)-0.01*min(x))
 
Z <- ZOSull(x,range.x,intKnots)
ZA <- (1-typeIsB)*Z
ZB <- typeIsB*Z

# Set up Z matrix for group specific curves:

numIntKnotsGrp <- 10
intKnotsGrp <- quantile(unique(x),
        seq(0,1,length = numIntKnotsGrp+2))[-c(1,numIntKnotsGrp+2)]

Zgrp <- ZOSull(x,range.x,intKnotsGrp)
ZgrpA <- (1-typeIsB)*Zgrp
ZgrpB <- typeIsB*Zgrp

# Specify hyperparameters:

sigmaBeta <- 1e5  ;  AU <- 1e5
AuGbl <- 1e5  ;      AuGrp <- 1e5

# Let up dimension variables:

ncXbase <- ncol(Xbase)
ncZ <- ncol(Z) ; ncZgrp <- ncol(Zgrp)

# Set up model in Stan:

bernGrpSpecContrModel <-
'data
{
   int<lower=1> numObs;             int<lower=1> numGrp;
   int<lower=1> idnum[numObs];      int<lower=1> ncXbase;
   int<lower=1> ncZ;                int<lower=1> ncZgrp;
   real<lower=0> sigmaBeta;         real<lower=0> AU;
   real<lower=0> AuGbl;             real<lower=0> AuGrp;
   int<lower=0,upper=1> y[numObs];  
   matrix[numObs,ncXbase] Xbase;    matrix[numObs,ncXbase] XA;
   matrix[numObs,ncXbase] XB;       real x[numObs];
   matrix[numObs,ncZ] ZA;           matrix[numObs,ncZ] ZB;
   matrix[numObs,ncZgrp] ZgrpA;     matrix[numObs,ncZgrp] ZgrpB;
}
transformed data
{
   vector[4] zeroVec;
   zeroVec = rep_vector(0,4); 
}
parameters
{
   vector[ncXbase] beta;         vector[ncXbase] betaBvsA;
   vector[ncZ] uGblA;            vector[ncZ] uGblB;
   vector[2*ncXbase] U[numGrp];  matrix[numGrp,ncZgrp] uGrpA;
   matrix[numGrp,ncZgrp] uGrpB;  cov_matrix[2*ncXbase] SigmaGrp;
   vector[2*ncXbase] a;          real<lower=0> siguGblA;
   real<lower=0> siguGblB;       real<lower=0> siguGrp;
}
transformed parameters
{
   vector[numObs] fmean;      vector[numObs] fullMean;
   fmean = Xbase*beta + XB*betaBvsA + ZA*uGblA + ZB*uGblB;   
   fullMean = fmean + to_vector(U[idnum,3]).*col(XA,1)
                    + to_vector(U[idnum,4]).*col(XA,2)
                    + to_vector(U[idnum,1]).*col(XB,1)
                    + to_vector(U[idnum,2]).*col(XB,2)
                         + rows_dot_product(uGrpA[idnum],ZgrpA)
                         + rows_dot_product(uGrpB[idnum],ZgrpB);
}
model
{
   vector[2*ncXbase] scaleSigmaGrp;

   y ~ bernoulli_logit(fullMean);

   U ~ multi_normal(zeroVec,SigmaGrp);

   uGblA ~ normal(0,siguGblA);  uGblB ~ normal(0,siguGblB);
   
   to_vector(uGrpA) ~ normal(0,siguGrp);
   to_vector(uGrpB) ~ normal(0,siguGrp);
      
   a ~ inv_gamma(0.5,pow(AU,-2)); 
   for (k in 1:(2*ncXbase)) scaleSigmaGrp[k] = 4/a[k];
   SigmaGrp ~ inv_wishart(5,diag_matrix(scaleSigmaGrp));

   beta ~ normal(0,sigmaBeta);   betaBvsA ~ normal(0,sigmaBeta);
   siguGblA ~ cauchy(0,AuGbl); siguGblB ~ cauchy(0,AuGbl);
   siguGrp ~ cauchy(0,AuGrp);
}'
	  
# Set up data list:

allData <- list(numObs = numObs,numGrp = numGrp,ncZ = ncZ,
                ncXbase = ncXbase,ncZgrp = ncZgrp,AU = AU,AuGbl = AuGbl,
                AuGrp = AuGrp,sigmaBeta = sigmaBeta,
                x = x,y = y,idnum = idnum,ZA = ZA,ZB = ZB,
     	        ZgrpA = ZgrpA,ZgrpB = ZgrpB,Xbase = Xbase,XA = XA,XB = XB)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = bernGrpSpecContrModel,data = allData,
                         iter = 1,chains = 1)

# Perform Markov chain Monte Carlo:

stanObj <-  stan(model_code = bernGrpSpecContrModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 1,
                 fit = stanCompilObj)
		   
betaMCMC <- NULL  ;  betaBvsAMCMC <- NULL   
for (i in 1:2)
{
   charVar1 <- paste("beta[",as.character(i),"]",sep = "") 
   betaMCMC <- rbind(betaMCMC,extract(stanObj,charVar1,permuted = FALSE))
   
   charVar2 <- paste("betaBvsA[",as.character(i),"]",sep = "") 
   betaBvsAMCMC <- rbind(betaBvsAMCMC,extract(stanObj,charVar2,permuted = FALSE))
}

uGblAMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("uGblA[",as.character(k),"]",sep = "") 
   uGblAMCMC <- rbind(uGblAMCMC,extract(stanObj,charVar,permuted = FALSE))
}
uGblBMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("uGblB[",as.character(k),"]",sep = "") 
   uGblBMCMC <- rbind(uGblBMCMC,extract(stanObj,charVar,permuted = FALSE))
}
siguGblAMCMC <- extract(stanObj,"siguGblA",permuted = FALSE)
siguGblBMCMC <- extract(stanObj,"siguGblB",permuted = FALSE)
siguGrpMCMC <- extract(stanObj,"siguGrp",permuted = FALSE)

UMCMC <- NULL
for (i in 1:numGrp)
   for (j in 1:4)
   {
      charVar <- paste("U[",as.character(i),",",as.character(j),"]",sep = "")
      UMCMC <- rbind(UMCMC,as.vector(extract(stanObj,charVar,permuted = FALSE)))
   }

Sigma11MCMC <- extract(stanObj,"SigmaGrp[1,1]",permuted = FALSE)
Sigma22MCMC <- extract(stanObj,"SigmaGrp[2,2]",permuted = FALSE)
Sigma33MCMC <- extract(stanObj,"SigmaGrp[3,3]",permuted = FALSE)
Sigma44MCMC <- extract(stanObj,"SigmaGrp[4,4]",permuted = FALSE)

uGrpAMCMC <- vector("list",numGrp); uGrpBMCMC <- vector("list",numGrp)
for (i in 1:numGrp)
   for (j in 1:ncZgrp)
   {
      charVar <- paste("uGrpA[",as.character(i),",",as.character(j),"]",sep = "")
      uGrpAMCMC[[i]] <- rbind(uGrpAMCMC[[i]],extract(stanObj,charVar,permuted = FALSE))
      charVar <- paste("uGrpB[",as.character(i),",",as.character(j),"]",sep = "")
      uGrpBMCMC[[i]] <- rbind(uGrpBMCMC[[i]],extract(stanObj,charVar,permuted = FALSE))
   }
    
# Do lattice-type graphic showing fitted group-specific curves:

AptCol <- "dodgerblue"   ;   AlnCol <- "blue"
BptCol <- "deeppink"     ;   BlnCol <- "maroon"
cexVal <- 0.2 ; lwdVal <- 1.5

ng <- 101
xg <- seq(min(x),max(x),length = ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,range.x,intKnots)
Zggrp <- ZOSull(xg,range.x,intKnotsGrp)

fAgMCMC <- Xg%*%betaMCMC + Zg%*%uGblAMCMC
fBgMCMC <- Xg%*%(betaMCMC + betaBvsAMCMC) + Zg%*%uGblBMCMC

gAgMCMC <- vector("list",numGrp)
gBgMCMC <- vector("list",numGrp)
currAMCMC <- vector("list",numGrp)
currBMCMC <- vector("list",numGrp) 
meanCrvAList <- vector("list",numGrp)
credLowAList <- vector("list",numGrp)
credUppAList <- vector("list",numGrp)
meanCrvBList <- vector("list",numGrp)
credLowBList <- vector("list",numGrp)
credUppBList <- vector("list",numGrp)

inds0A <- (1:length(x))[(y == 0)&(typeIsB == 0)]
inds0B <- (1:length(x))[(y == 0)&(typeIsB == 1)]
inds1A <- (1:length(x))[(y == 1)&(typeIsB == 0)]
inds1B <- (1:length(x))[(y == 1)&(typeIsB == 1)]

inds0Alist <- vector("list",numGrp)
inds0Blist <- vector("list",numGrp)
inds1Alist <- vector("list",numGrp)
inds1Blist <- vector("list",numGrp)

nVec <- as.vector(table(idnum))
reBlockInds <- vector("list",length = numGrp)
currStt <- 1

for (i in 1:numGrp)
{
   currEnd <- currStt + nVec[i] - 1
   reBlockInds[i] <- list(currStt:currEnd)
   currStt <- currEnd + 1

   currPanelInds <- (1:length(x))[idnum == i]
   inds0Alist[[i]] <- intersect(currPanelInds,inds0A) 
   inds0Blist[[i]] <- intersect(currPanelInds,inds0B)
   inds1Alist[[i]] <- intersect(currPanelInds,inds1A) 
   inds1Blist[[i]] <- intersect(currPanelInds,inds1B)
   
   gAgMCMC[[i]] <- (Xg%*%rbind(UMCMC[(4*i-1),],UMCMC[(4*i),])  
                    + Zggrp%*%uGrpAMCMC[[i]])
   gBgMCMC[[i]] <- (Xg%*%rbind(UMCMC[(4*i-3),],UMCMC[(4*i-2),])  
                    + Zggrp%*%uGrpBMCMC[[i]])
   currAMCMC[[i]] <- fAgMCMC + gAgMCMC[[i]] 
   meanCrvAList[[i]] <- apply(plogis(currAMCMC[[i]]),1,mean)
   credLowAList[[i]] <- apply(plogis(currAMCMC[[i]]),1,quantile,0.025)
   credUppAList[[i]] <- apply(plogis(currAMCMC[[i]]),1,quantile,0.975)

   currBMCMC[[i]] <- fBgMCMC + gBgMCMC[[i]]    
   meanCrvBList[[i]] <- apply(plogis(currBMCMC[[i]]),1,mean)
   credLowBList[[i]] <- apply(plogis(currBMCMC[[i]]),1,quantile,0.025)
   credUppBList[[i]] <- apply(plogis(currBMCMC[[i]]),1,quantile,0.975)
}

allDataDF <- data.frame(x,y,idnum,typeIsB)
xOuter <- x; yOuter <- y

figFit <- xyplot(y~x|idnum,data = allDataDF,
                 layout = c(5,5),main = "",
                 xlab = list("log(initial size (cm) +1)",cex = 1.5),
                 ylab = list("indicator that coral is alive",cex = 1.5),
                 cex.lab = 4,scales = list(cex = 1.25),
                 strip = FALSE,xlim = c(0,4),ylim = c(-0.3,1.3),
                 key  =  list(title = "",
                      columns  =  2, 
                      points  =  list(pch  =  rep(16,2),
                      col  =  c(AptCol,BptCol),cex = 4*cexVal),
                      text  =  list(c("Pocillopora",
                                    "Porites"),cex = 1.8)),
                 panel = function(x,y,subscripts,groups)
                 {
                    panel.grid()
                    iGrp <- panel.number()

                    xCurr0A <- xOuter[inds0Alist[[iGrp]]]
                    yCurr0A <- runif(length(xCurr0A),-0.3,0)

                    xCurr1A <- xOuter[inds1Alist[[iGrp]]]
                    yCurr1A <- runif(length(xCurr1A),1,1.3)

                    xCurr0B <- xOuter[inds0Blist[[iGrp]]]
                    yCurr0B <- runif(length(xCurr0B),-0.3,0)

                    xCurr1B <- xOuter[inds1Blist[[iGrp]]]
                    yCurr1B <- runif(length(xCurr1B),1,1.3)

                    panel.xyplot(xCurr0A,yCurr0A,col = AptCol,cex = cexVal)
                    panel.xyplot(xCurr1A,yCurr1A,col = AptCol,cex = cexVal)
                    panel.xyplot(xCurr0B,yCurr0B,col = BptCol,cex = cexVal)
                    panel.xyplot(xCurr1B,yCurr1B,col = BptCol,cex = cexVal)

                    panel.xyplot(xg,meanCrvAList[[iGrp]],col = AlnCol,
                                 type = "l",lwd = lwdVal)
                    panel.xyplot(xg,credLowAList[[iGrp]],col = AlnCol,
                                 type = "l",lty = 2,lwd = lwdVal)
                    panel.xyplot(xg,credUppAList[[iGrp]],col = AlnCol,
                                 type = "l",lty = 2,lwd = lwdVal)
                    panel.xyplot(xg,meanCrvBList[[iGrp]],col = BlnCol,
                                 type = "l",lwd = lwdVal)
                    panel.xyplot(xg,credLowBList[[iGrp]],col = BlnCol,
                                 type = "l",lty = 2,lwd = lwdVal)
                    panel.xyplot(xg,credUppBList[[iGrp]],col = BlnCol,
                                 type = "l",lty = 2,lwd = lwdVal)
                    panel.xyplot(xg,rep(0,length(xg)),col = "navy",
                                 type = "l",lwd = lwdVal,lty = 3)
                    panel.xyplot(xg,rep(1,length(xg)),col = "navy",
                                 type = "l",lwd = lwdVal,lty = 3)

                 })
print(figFit)

# Do MCMC summary plot to check convergence:

ContgMCMC <- fBgMCMC - fAgMCMC

indQ1 <- length(xg[xg<quantile(x,0.25)])
indQ2 <- length(xg[xg<quantile(x,0.50)])
indQ3 <- length(xg[xg<quantile(x,0.75)])
ContQ1 <- ContgMCMC[indQ1,]
ContQ2 <- ContgMCMC[indQ2,]
ContQ3 <- ContgMCMC[indQ3,]

MCMClist <- list(cbind(ContQ1,ContQ2,ContQ3))
parNamesVal <- list(c("contrast funct.","at 1st quartile","of log(init. size +1)"),
                    c("contrast funct.","at 2nd quartile","of log(init. size +1)"),
                    c("contrast funct.","at 3rd quartile","of log(init. size +1)"))

summMCMC(MCMClist,parNames = parNamesVal)

# Do contrast curve plot:

ContgMean <- apply(ContgMCMC,1,mean)
ContgLower <- apply(ContgMCMC,1,quantile,0.025)
ContgUpper <- apply(ContgMCMC,1,quantile,0.975)

cex.labVal <- 1.5 ; cex.axisVal <- 1.5
ylimVal <- c(min(c(ContgLower,ContgUpper)),3)
par(mai = c(1,0.95,0.1,0.1))
plot(xg,ContgMean,bty = "l",type = "n",xlab = "log(initial size (cm) + 1)",
     ylab = "log(odds ratio)",ylim = ylimVal,cex.lab = cex.labVal,cex.axis = cex.axisVal)

polygon(c(xg,rev(xg)),c(ContgLower,rev(ContgUpper)),
         col = "palegreen",border = FALSE)
lines(xg,ContgMean,col = "darkgreen",lwd = 2)
abline(0,0,col = "slateblue")
nSub <- 1000
subInds <- sample(1:length(x),nSub,replace = FALSE)
rug(x[subInds],col = "dodgerblue")

############ End of coralAna ###########








