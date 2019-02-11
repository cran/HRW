########## R script: returnsDrvFig ##########

# For creation of the figure from Section 5.4 
# figure comparing the slope of the fitted
# surface of the varying-coefficient model
# with that of the bivariate nonparametric 
# regression model for data on Standard & Poor's
# 500 excess log-returns.

# Last changed: 08 NOV 2017

# Load required packages:

library(mgcv) ; library(HRW)

# Obtain the excess log-returns data:

data(capm) ; n <- dim(capm)[1]
riskfree <- capm$Close.tbill[2:n]/(100*365)
elrGE <- diff(log(capm$Close.ge)) - riskfree
elrSP500 <- diff(log(capm$Close.sp500))  - riskfree

# Create the time variable `t':

dayNums <- (1:(n-1))/(n-1) 
startTime <- 1993 + 11/12 ; endTime <- 2003 + 3/12
t <- startTime + (endTime - startTime)*dayNums

# Obtain the varying-coefficient model and bivariate
# nonparametric regression fits:

fitVCM <- gam(elrGE ~ s(t) + s(t,by = elrSP500),
              method = "REML")

fitBivNPR <- gam(elrGE ~ te(t,elrSP500,k = rep(25,2)),
                 method = "REML")

# Plot the slope component of the varying-coefficient model fit:

par(mfrow = c(1,2),mai = c(0.9,0.9,0.8,0.05))
ylimVal <- c(0.8,1.5)
plot(fitVCM,xlab = "time (year)",ylab = "slope",rug = FALSE, 
     main = "varying-coefficient model",shade = TRUE, 
     shade.col = "palegreen",select = 2,
     col = "darkgreen",bty = "l",lwd = 2,
     col.main = "navy",cex.main = 0.85,ylim = ylimVal)
rug(t,col = "dodgerblue")

# Obtain and plot numerical derivative of the bivariate nonparametric
# regression estimate with respect to `t':

ngt <- 101 ; tgrid <- seq(min(t),max(t),length = ngt)
elrSP500Val <- mean(elrSP500) ; hVal <- 0.000001
fMinushGrid <- predict(fitBivNPR,data.frame(elrSP500 
                       = rep((0 - hVal),ngt),t = tgrid))
fPlushGrid <- predict(fitBivNPR,data.frame(elrSP500 
                       = rep((0 + hVal),ngt),t = tgrid))
slopeGrid <- (fPlushGrid - fMinushGrid)/(2*hVal)
plot(tgrid,slopeGrid,bty = "l",col = "darkgreen", 
     lwd = 2,type = "l",xlab = "time (year)",ylab = "slope",
     main = "bivariate nonparametric regression model",
     col.main = "navy",cex.main = 0.85,ylim = ylimVal)
rug(t,col = "dodgerblue")

############ End of returnsDrvFig ############
