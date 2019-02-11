########## R script: absorbFuncDrvs ##########

# For plotting absorbance functions and their
# first two derivatives for the `Fat Content of
# Meat Samples' example in Section 6.3.2.

# Last changed: 25 OCT 2017

# Read in functional data:

library(fda) ; library(fda.usc) ; data(tecator)
X <- tecator$absorp.fdata$data

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

# Set plotting layout:

par(mfrow = c(1,3),mai = c(0.6,0.6,0.2,0.15))

# Plot the first ten absorbance functions:

plot(wavelength,X[1,],type = "l",ylim = range(X[1:10,]), 
     xlab = "wavelength (nanometers)",
     ylab= "absorbance",bty = "l",
     main = "absorbance functions",
     cex.lab = 1.4,cex.axis = 1.5, 
     col = "darkgreen",cex.main = 1.5,
     col.main = "navy")
for (i in 2:10) lines(wavelength,X[i,],col = i)

# Plot their first derivatives:

plot(wavelength,Xderiv[1,],type = "l",
     ylim = range(Xderiv[1:10,]),
     xlab = "wavelength (nanometers)",
     ylab = "1st derivative of absorbance",
     main = "first derivatives",
     bty = "l",cex.lab = 1.4,cex.axis = 1.5,
     col = "darkgreen",cex.main = 1.5,col.main = "navy")
for (i in 2:10) lines(wavelength,Xderiv[i,],col = i)
abline(h = 0,col = "slateblue")

# Plot their second derivatives:

plot(wavelength,Xderiv2[1,],type = "l",
     ylim = range(Xderiv2[1:10,]),
     xlab = "wavelength (nanometers)",
     ylab = "2nd derivative of absorbance",
     main = "second derivatives",
     bty = "l",cex.lab = 1.4,cex.axis = 1.5,
     col = "darkgreen",cex.main = 1.5,col.main = "navy")
for (i in 2:10) lines(wavelength,Xderiv2[i,],col = i)
abline(h = 0,col = "slateblue")

############ End of absorbFuncDrvs ############
