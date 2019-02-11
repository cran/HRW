########## R script: absorbBivarFigs ##########

# For obtaining the bivariate figure for the 
# fat content of meat samples absorbance functions 
# example in Section 6.4.1.

# Last changed: 08 NOV 2017

# Read in functional data corresponding absorbance of infrared light
# at 100 equally spaced wavelengths from 850 to 1,050 nanometers
# for 215 finely chopped meat samples and response data corresponding
# to the percentage of fat:

library(refund) ; library(fda) 
library(fda.usc) ; data(tecator)
y <- tecator$y$Fat ; X <- tecator$absorp.fdata$data

# Create array containing wavelength values:

wavelength <- seq(850,1050,length = 100)

# Obtain first and second derivatives of the absorbance functions:

bbt <- create.bspline.basis(range(wavelength), nbasis=20)  
Xfd <- smooth.basisPar(wavelength,t(X),bbt,2,1e-9)
Xderiv <- t(eval.fd(wavelength,Xfd$fd,1))  
Xderiv2 <- t(eval.fd(wavelength,Xfd$fd,2)) 

# Obtain the scalar-on-function additive model fit to the raw data:

fitPFRnoTrans <- pfr(y ~ af(Xderiv2,argvals=wavelength, 
   bs = "ps",k = c(7,7),m = list(c(2,2),c(2,2)),
   Qtransform = FALSE), method = "REML")

# Obtain the scalar-on-function additive model fit to the
# empirical cumulative distribution function transformed data:

fitPFRtrans <- pfr(y ~ af(Xderiv2,argvals = wavelength,
   bs = "ps",k = c(7,7),m = list(c(2,2),c(2,2)),
   Qtransform = TRUE),method = "REML")

# Set plotting layout:

par(mfrow=c(2,2),mai=c(0.8,0.7,0.3,0.15))

# Plot the raw data fit with `rug = FALSE':

plot(fitPFRnoTrans,scheme = 2,
     main = "Qtransform = FALSE,rug = FALSE", 
     xlab = "t",ylab = "x",bty = "l",cex.lab = 1.5, 
     cex.axis = 1.5,cex.main = 1.6,col.main = "navy",
     rug = FALSE)

# Plot the raw data fit with `rug = TRUE':

plot(fitPFRnoTrans,scheme = 2,
     main = "Qtransform = FALSE,rug = TRUE",
     xlab = "t",ylab = "x",bty = "l",
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.6,
     col.main = "navy",rug = TRUE)

# Plot the transformed data fit with `rug = FALSE':

plot(fitPFRtrans,scheme = 2, 
     main = "Qtransform = TRUE,rug = FALSE",
     xlab = "t",ylab = "x",bty = "l",
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.6,
     col.main = "navy",rug = FALSE,Qtransform = TRUE)

# Plot the transformed data fit with `rug = TRUE':

plot(fitPFRtrans,scheme = 2,
     main = "Qtransform = TRUE,rug = TRUE",
     xlab = "t",ylab = "x",bty = "l",
     cex.lab = 1.5,cex.axis = 1.5,cex.main = 1.6,
     col.main = "navy",rug = TRUE,Qtransform = TRUE)

############ End of absorbBivarFigs ############
