########## R script: MichIncLMS ##########

# For doing LMS quantile regression to two variables
# from the 1987 cross-section of the Michigan Panel
# Study of Income Dynamics.

# Last changed: 19 JUN 2017

# Load required package:

library(VGAM)

# Read in data:

library(Ecdat)
data(Workinghours)
x <- Workinghours$age 
y <- Workinghours$income/10

# Set up plotting infrastructure:

cexVal <- 0.4 ; cex.labVal <- 2.6 ; lwdVal <- 2
cex.axisVal <- 2.2 ; cex.mainVal <- 2.5
col.mainVal <- "navy"
pointCol <- "dodgerblue" ; lineCol <- c("darkgreen","indianred3")
ylabString <- "other household income ('000 $US)"
mainStrings <- c("full data view","zoomed view")
ng <- 201
xg <- seq(min(x),max(x),length = ng)
par(mfrow = c(1,2),mai = c(0.9,0.9,0.54,0.05)) 

# Set quantile vector and smoothing parameters:

quantVec <- c(1,5,25,50,75,95,99)/100
colVec <- c("darkmagenta","navy","blue","darkgreen",
             "gold","darkorange","red")

# Remove the tiny proportion of observations that are negative
# and also work with scaled version of response data (divided by 100)
# for stability.

npi <- which(y <= 0)
yLMS <- y[-npi]/100 ; xLMS <- x[-npi]

fitLMSinit <- vgam(yLMS ~ s(xLMS,df = 4),lms.bcn(zero = c(1,3)),maxit = 50)

fitLMS  <- vgam(yLMS ~ s(xLMS,df = c(3,10,3)),lms.bcn(zero = NULL), 
                etastart  =  predict(fitLMSinit),maxit = 500)

xgDF <- data.frame(xLMS = xg)
LMSgrids <- predict(fitLMS,xgDF)

lambdag <- LMSgrids[,1]
mug <- LMSgrids[,2]
logSigmag <- LMSgrids[,3]

fitLMSg <- vector("list",length(quantVec))
for (iq in 1:length(quantVec))
   fitLMSg[[iq]] <- mug*(1 + lambdag*exp(logSigmag)*qnorm(quantVec[iq]))^(1/lambdag) 

par(mfrow = c(1,2),mai = c(0.9,0.9,0.54,0.05)) 

# First do plots with actual data:

plot(x,y,cex = cexVal,bty = "l",xlab = "wife's age in years",
     ylab = ylabString,
     cex.lab = cex.labVal,col = pointCol,main = mainStrings[1],
     cex.main = cex.mainVal,cex.axis = cex.axisVal,col.main = col.mainVal)

legend("topleft",legend = c("99% quantile","95% quantile","75% quantile",
                          "50% quantile","25% quantile","5% quantile",
                          "1% quantile"),lwd = rep(lwdVal,7),col = rev(colVec),
                          cex = 1.5)

for (iq in 1:length(quantVec))
   lines(xg,100*fitLMSg[[iq]],col = colVec[iq],lwd = lwdVal)

plot(x,y,cex = cexVal,bty = "l",xlab = "wife's age in years",
     ylab = ylabString,
     ylim = c(0,200),cex.lab = cex.labVal,col = pointCol,main = mainStrings[2],
     cex.main = cex.mainVal,cex.axis = cex.axisVal,col.main = col.mainVal)

for (iq in 1:length(quantVec))
   lines(xg,100*fitLMSg[[iq]],col = colVec[iq],lwd = lwdVal)

# Now do (lambda,mu,sigma) plots:

par(mfrow = c(1,3))
plot(xg,lambdag,bty = "l",col = "darkgreen",type = "l",main = "estimated lambda function")
plot(xg,mug,bty = "l",col = "darkgreen",type = "l",main = "estimated mu function")
plot(xg,exp(logSigmag),bty = "l",col = "darkgreen",type = "l",main = "estimated sigma function")

############ End of MichIncLMS ############

