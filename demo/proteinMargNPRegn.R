########## R script: proteinMargNPRegn ##########

# For conducting a marginal nonparametric regression 
# for the female-only protein data:

# Last changed: 22 AUG 2017

# Set flag for code compilation (needed if 
# running script first time in current session):

compileCode <- TRUE

# Set MCMC sample size paramaters:

nWarm <- 1000
nKept <- 1000
nThin <- 1

# Load required packages:

library(HRW)   ;   library(rstan)   ;   library(lattice)

# Extract regression data: 

data(protein)
femInds <- (1:nrow(protein))[protein$female==1]
femProtein <- protein[femInds,]
xOrig <- femProtein$BMI
yOrig <- femProtein$proteinBioM
m <- 130 ; n <- 2

# Work with standardized data for Bayesian analysis:

mean.x <- mean(xOrig)   ;   sd.x <- sd(xOrig)
mean.y <- mean(yOrig)   ;   sd.y <- sd(yOrig)
x <- (xOrig - mean.x)/sd.x
y <- (yOrig - mean.y)/sd.y

# Convert x and y to matrix format:

xMat <- matrix(x,m,n,byrow=TRUE)
yMat <- matrix(y,m,n,byrow=TRUE)

# Set up basic variables for the spline component:

range.x.Val <- c(1.05*min(x)-0.05*max(x),1.05*max(x)-0.05*min(x))
numIntKnots <- 15
intKnots <- seq(range.x.Val[1],range.x.Val[2],
                length=(numIntKnots+2))[-c(1,numIntKnots+2)]

# Obtain the spline component of the Z matrix:

Z <- ZOSull(x,intKnots=intKnots,range.x=range.x.Val)
ncZ <- ncol(Z)

# Set hyperparameters:

sigmaBeta <- 1e5 ; Au <- 1e5 ; ASigma <- 1e5

# Specify model in Stan:

margNPRegnModel <- 
'data
{
   int<lower=1> m;          int<lower=1> n;
   int<lower=1> ncZ;        real<lower=0> sigmaBeta;    
   real<lower=0> Au;        real<lower=0> ASigma;
   matrix[m,n] xMat;        matrix[m,n] yMat;
   matrix[m*n,ncZ] Z;
}
parameters
{
   real beta0;              real beta1;
   vector[ncZ] u;           real<lower=0> sigma;    
   cov_matrix[n] Sigma;     vector[n] a;
}
transformed parameters
{
   matrix[m,n] meanFuncPart1;
   matrix[m,n] meanFuncPart2;
   matrix[m,n] meanFuncFull;

   meanFuncPart1 = beta0 + beta1*xMat;

   for (i in 1:m) 
      for (j in 1:n) 
         meanFuncPart2[i,j] = dot_product(u,Z[(i-1)*n+j]);

   meanFuncFull = meanFuncPart1 + meanFuncPart2 ;  
}
model
{
   vector[n] scaleSigma;

   for (i in 1:m) 
      yMat[i] ~ multi_normal(meanFuncFull[i],Sigma);	
   	
   u ~ normal(0,sigma);
   
   a ~ inv_gamma(0.5,pow(ASigma,-2));     
   for (j in 1:n) scaleSigma[j] = 4/a[j];
   Sigma ~ inv_wishart(n+1,diag_matrix(scaleSigma));

   beta0  ~ normal(0,sigmaBeta) ; beta1 ~ normal(0,sigmaBeta);
   sigma ~ cauchy(0,Au);
}'

# Set up data list:

allData <- list(m = m,n = n,ncZ = ncZ,xMat = xMat,yMat = yMat,Z = Z,
                sigmaBeta = sigmaBeta,Au = Au,ASigma = ASigma)

# Compile code for model if required:

if (compileCode)
   stanCompilObj <- stan(model_code = margNPRegnModel,data = allData,
                         iter = 1,chains = 1)

# Perform MCMC:

stanObj <-  stan(model_code = margNPRegnModel,data = allData,warmup = nWarm,
                 iter = (nWarm + nKept),chains = 1,thin = nThin,refresh = 100,
                 fit = stanCompilObj)

# Save and extract relevant MCMC samples:

beta0MCMC <- as.vector(extract(stanObj,"beta0",permuted = FALSE))
beta1MCMC <- as.vector(extract(stanObj,"beta1",permuted = FALSE))
betaMCMC <- rbind(beta0MCMC,beta1MCMC)

uMCMC <- NULL
for (k in 1:ncZ)
{
   charVar <- paste("u[",as.character(k),"]",sep="") 
   uMCMC <- rbind(uMCMC,extract(stanObj,charVar,permuted=FALSE))
}

sigmaMCMC <- as.vector(extract(stanObj,"sigma",permuted=FALSE))
Sigma11MCMC <- extract(stanObj,"Sigma[1,1]",permuted = FALSE)
Sigma12MCMC <- extract(stanObj,"Sigma[1,2]",permuted = FALSE)
Sigma22MCMC <- extract(stanObj,"Sigma[2,2]",permuted = FALSE)

sigmaMCMCOrig <- sigmaMCMC*(sd.y^2)
Sigma11MCMCOrig <- Sigma11MCMC*(sd.y^2)
Sigma12MCMCOrig <- Sigma12MCMC*(sd.y^2)
Sigma22MCMCOrig <- Sigma22MCMC*(sd.y^2)

# Compute the curve estimate:

ng <- 201
xg <- seq(range.x.Val[1],range.x.Val[2],length=ng)
Xg <- cbind(rep(1,ng),xg)
Zg <- ZOSull(xg,intKnots=intKnots,range.x=range.x.Val) 
fhatMCMC <- Xg%*%betaMCMC + Zg%*%uMCMC

# Convert fhatMCMC matrix to original scale:

fhatMCMCorig <- mean.y + sd.y*fhatMCMC

# Obtain plotting grids on the orginal scale:

xgOrig <- mean.x + sd.x*xg
fhatgOrig <- apply(fhatMCMCorig,1,mean)
credLowergOrig <- apply(fhatMCMCorig,1,quantile,0.025)
credUppergOrig <- apply(fhatMCMCorig,1,quantile,0.975)

# Obtain vertical slices at the quartiles of the predictor:

indQ1 <- length(xgOrig[xgOrig<quantile(xOrig,0.25)])
indQ2 <- length(xgOrig[xgOrig<quantile(xOrig,0.50)])
indQ3 <- length(xgOrig[xgOrig<quantile(xOrig,0.75)])
fhatOrigQ1 <- fhatMCMCorig[indQ1,]
fhatOrigQ2 <- fhatMCMCorig[indQ2,]
fhatOrigQ3 <- fhatMCMCorig[indQ3,]

# Make an an MCMC summary plot:

parmsMCMC <- list(cbind(fhatOrigQ1,fhatOrigQ2,fhatOrigQ3,
                        Sigma11MCMCOrig,Sigma12MCMCOrig,Sigma22MCMCOrig))
parNamesVec <- list(c("mean function","at 1st quartile","of BMI"),
                    c("mean function","at 2nd quartile","of BMI"),
                    c("mean function","at 3rd quartile","of BMI"),
                    expression(Sigma[11]),expression(Sigma[12]),
                    expression(Sigma[22]))
summMCMC(parmsMCMC,parNames=parNamesVec,numerSummCex=1.3)

# Do fit plot, starting with setting up objects for an L-shaped axis
# in xyplot():
  
stripSuppression <- list()
stripSuppression$axis.line$col <- NA
stripSuppression$strip.border$col <- NA
stripSuppression$strip.background$col <- NA

axisL <- function(side,...,line.col)
{
   if (side %in% c("bottom", "left")) 
   {
      col <- trellis.par.get("axis.text")$col
      axis.default(side,...,line.col=col)
      if (side=="bottom")
         grid::grid.lines(y=0)
      if (side=="left")
         grid::grid.lines(x=0)
   }
}

cex.labVal <- 1.8  ; cex.axisVal <- 1.6
mNPrFit <- xyplot(proteinBioM~BMI,groups=idnum,data=femProtein,type="b",
                  subscripts=TRUE,axis=axisL,par.settings=stripSuppression,
                  xlab=list(label="body mass index",cex=cex.labVal),
                  ylab=list(label="logarithm of protein biomarker",cex=cex.labVal),
                  scales=list(x=list(cex=cex.axisVal),y=list(cex=cex.axisVal)),
                  panel=function(x,y,subscripts,groups)
                  {
                     panel.polygon(c(xgOrig,rev(xgOrig)),
                                   c(credLowergOrig,rev(credUppergOrig)),
                                   col="palegreen",border=FALSE)
                     panel.xyplot(xgOrig,fhatgOrig,lwd=2,type="l",col="darkgreen")
                     panel.superpose(x,y,subscripts,groups,type="b")
                  }
                  )
plot(mNPrFit)

########## End of proteinMargNPRegn ##########


