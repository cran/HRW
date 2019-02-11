########## R script: MultIncMultQSS ##########

# For doing multiple quantile smoothing spline fits
# to two variables from the 1987 cross-section
# of the Michigan Panel Study of Income Dynamics.

# Last changed: 20 APR 2017

# Load required package:

library(quantreg)

# Read in data:

library(Ecdat)
data(Workinghours)
x <- Workinghours$age 
y <- Workinghours$income/10

# Set up plotting infrastructure:

cexVal <- 0.4 ; cex.labVal <- 1 ; lwdVal <- 2
cex.axisVal <- 1 ; cex.mainVal <- 1
col.mainVal <- "navy"
pointCol <- "dodgerblue" ; lineCol <- c("darkgreen","indianred3")
ylabString <- "other household income ('000 $US)"
mainStrings <- c("full data view","zoomed view")
par(mfrow = c(1,2),mai = c(1,0.8,0.3,0.05))
ng <- 201
xg <- seq(min(x),max(x),length = ng)

# Set quantile vector and smoothing parameters:

quantVec <- c(0.01,0.05,0.25,0.5,0.75,0.95,0.99)
colVec <- c("darkmagenta","navy","blue","darkgreen",
                "gold","darkorange","red")
lambdaVal <- 3.5

fitQSS <- vector("list",length(quantVec))
fitQSSg <- vector("list",length(quantVec))
for (iq in 1:length(quantVec))
{
   fitQSS[[iq]] <- rqss(y ~ qss(x,lambda = lambdaVal),tau = quantVec[iq])
   xgDF <- data.frame(x = xg)
   fitQSSg[[iq]] <- predict(fitQSS[[iq]],xgDF)
}

# First do plots with actual data:

plot(x,y,cex = cexVal,bty = "l",xlab = "wife's age in years",
     ylab = ylabString,cex.lab = cex.labVal,col = pointCol,main = mainStrings[1],
     cex.main = cex.mainVal,cex.axis = cex.axisVal,col.main = col.mainVal)

legend("topleft",legend = c("99% quantile","95% quantile","75% quantile",
                          "50% quantile","25% quantile","5% quantile",
                          "1% quantile"),lwd = rep(lwdVal,7),col = rev(colVec))

for (iq in 1:length(quantVec))
   lines(xg,fitQSSg[[iq]],col = colVec[iq],lwd = lwdVal)

# Now do plot with zoomed data:

plot(x,y,cex = cexVal,bty = "l",xlab = "wife's age in years",
     ylab = ylabString,ylim = c(0,200),cex.lab = cex.labVal,col = pointCol,
     main = mainStrings[2],cex.main = cex.mainVal,cex.axis = cex.axisVal,
     col.main = col.mainVal)

for (iq in 1:length(quantVec))
     lines(xg,fitQSSg[[iq]],col = colVec[iq],lwd = lwdVal)

############ End of MultIncMultQSS ############

