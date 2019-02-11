########## R script: absorbScalOnFuncAna ##########

# For illustration of scalar-on-function regression
# using data on meat content and absorbance spectra.

# Last changed: 06 NOV 2017

# Load in required packages:

library(fda) ; library(fda.usc) ; library(refund)

# Set up plotting layout:

par(mfrow = c(1,3),mai = c(0.6,0.9,0.30,0.1))

# Load in scalar response (y = percentages of fat, water and
# protein) and functional predictor (X = absorbance functions) 
# data:

data(tecator)
y <- tecator$y$Fat ; X <- tecator$absorp.fdata$data

# Create an array containing the wavelength grid
# corresponding to the absorbance function ordinates
# in X:

wavelength <- seq(850,1050,length = 100)

# Create the matrices `Xderiv' and `Xderiv2' containing 
# the first and second derivatives of the absorbance
# functions:

bbt <- create.bspline.basis(rangeval = range(wavelength),
                            nbasis = 20)   
Xfd <- smooth.basisPar(wavelength,t(X),bbt,2,1e-9)  
Xderiv <- t(eval.fd(wavelength,Xfd$fd,1))  
Xderiv2 <- t(eval.fd(wavelength,Xfd$fd,2)) 

# Perform the scalar-on-function regression for the
# absorbance functions using refund:::pfr():

fitlf0 <- pfr(y ~ lf(X,argvals = wavelength),method = "REML")

# Plot the estimated coefficient coefficient function:

plot(fitlf0,shade = TRUE,col = "darkgreen",
     shade.col = "palegreen",xlab = "wavelength (nanometers)",
     ylab = expression(widehat(beta)[1](t)),bty = "l",
     main = "zeroth derivative predictors",cex.lab = 1.5,
     cex.axis = 1.5,cex.main = 1.5,col.main = "navy")
abline(h = 0,col = "slateblue")

# Perform the scalar-on-function regression for the
# first derivative of the absorbance functions:

fitlf1 <- pfr(y ~ lf(Xderiv,argvals = wavelength),
              method = "REML")

# Plot the fit:

plot(fitlf1,shade = TRUE,col = "darkgreen",
     shade.col = "palegreen",xlab = "wavelength (nanometers)",
     ylab = expression(widehat(beta)[1](t)), 
     bty = "l",main = "first derivative predictors",
     cex.lab = 1.5,cex.axis = 1.5,
     cex.main = 1.5,col.main = "navy")
abline(h = 0,col = "slateblue")

# Perform the scalar-on-function regression for the
# second derivative of the absorbance functions:

fitlf2 <- pfr(y ~ lf(Xderiv2,argvals = wavelength),
              method = "REML")

# Plot the fit:

plot(fitlf2,shade = TRUE,col = "darkgreen",
     shade.col = "palegreen",xlab = "wavelength (nanometers)",
     ylab = expression(widehat(beta)[1](t)), 
     bty = "l",main = "second derivative predictors",
     cex.lab = 1.5,cex.axis = 1.5,
     cex.main = 1.5,col.main = "navy")
abline(h = 0,col = "slateblue")

############ End of absorbScalOnFuncAna ############
