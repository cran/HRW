########## R script: OFPGAMfit ##########

# For producing plot of the estimated smooth function
# components for the Poisson generalized additive model
# fit to the physician office visits data on both the
# link scale and the response scale.

# Last changed: 14 MAR 2018

# Load required packages:

library(HRW) ; library(mgcv) ; library(Ecdat) 

# Obtain data to be used for generalized additive 
# model analysis:

data(OFP) ; OFPforAna <- OFP
OFPforAna$age <- 10*OFPforAna$age
OFPforAna <- OFPforAna[OFPforAna$age <= 95,]

# Obtain generalized additive model fit:

fitGAMOFP <- gam(ofp ~ s(age) + s(school) + adldiff +  black 
                 + sex + maried + privins + medicaid + region + hlth,
	         family = poisson, scale = -1, data = OFPforAna)

# Plot the estimated smooth function components fits on both 
# the link scale and the response scale:

par(mfrow = c(2,2),lwd = 2,bty = "l",cex.lab = 1.5,
    cex.main = 1.5,col.main = "navy",cex.axis = 1.5,
    mai = c(0.9,1.2,0.35,0.05))

plot(fitGAMOFP,select = 1,shade = TRUE,shade.col = "palegreen",
     main = "link scale",xlab = "age in years",ylab="",rug=FALSE)
rug(jitter(OFPforAna$age),col = "dodgerblue")
mtext("effect on the log of mean no.",line = 4,side = 2,cex = 1.25)
mtext("of physician office visits",line = 3,side = 2,cex = 1.25)

plot(fitGAMOFP,select = 2,shade = TRUE,shade.col = "palegreen",
     main = "link scale",xlab = "number of years of education",
     ylab = "",rug = FALSE)
rug(jitter(OFPforAna$school),col = "dodgerblue")
mtext("effect on the log of mean no.",line = 4,side = 2,cex = 1.25)
mtext("of physician office visits",line = 3,side = 2,cex = 1.25)

# Set function for obtaining the mode of a vector:

modalValue <- function(x)
   return(unique(x)[which.max(tabulate(match(x,unique(x))))])

# Set grids for `age' and `school':

ng <- 401 ; ageg <- seq(66,95,length = ng) ; schoolg <- seq(0,18,length = ng)

# Obtain and plot slice of the mean surface
# in the `age' direction corresponding to the modal
# values of categorical predictors and the mean of 
# other continuous predictors:

newdataDF <- data.frame(age = ageg,
                        school = mean(OFPforAna$school),
                        adldiff = modalValue(OFPforAna$adldiff),
                        black  = modalValue(OFPforAna$black),
                        sex  = modalValue(OFPforAna$sex),
                        maried  = modalValue(OFPforAna$maried),
                        privins  = modalValue(OFPforAna$privins),
                        medicaid  = modalValue(OFPforAna$medicaid),
                        region  = modalValue(OFPforAna$region),
                        hlth  = modalValue(OFPforAna$hlth))

predObjage <- predict(fitGAMOFP,newdata = newdataDF,
                      type = "response",se.fit = TRUE)
etahatageg <- predObjage$fit
lowageg <- etahatageg - qnorm(0.975)*predObjage$se.fit
uppageg <- etahatageg + qnorm(0.975)*predObjage$se.fit

plot(0,type = "n",xlim = range(OFPforAna$age),ylim = c(2.5,8),xlab  =  "age in years",
     ylab  =  "mean no. of physician office visits",main  =  "response scale")
polygon(c(ageg,rev(ageg)),c(lowageg,rev(uppageg)),col = "palegreen",border = FALSE)
lines(ageg,etahatageg,col = "darkgreen",lwd = 2)
rug(jitter(OFPforAna$age),col = "dodgerblue",quiet = TRUE)

# Obtain and plot slice of the mean surface
# in the `school' direction corresponding to the modal
# values of categorical predictors and the mean of 
# other continuous predictors:

newdataDF <- data.frame(age = mean(OFPforAna$age),
                        school = schoolg,
                        adldiff = modalValue(OFPforAna$adldiff),
                        black  = modalValue(OFPforAna$black),
                        sex  = modalValue(OFPforAna$sex),
                        maried  = modalValue(OFPforAna$maried),
                        privins  = modalValue(OFPforAna$privins),
                        medicaid  = modalValue(OFPforAna$medicaid),
                        region  = modalValue(OFPforAna$region),
                        hlth  = modalValue(OFPforAna$hlth))

predObjschool <- predict(fitGAMOFP,newdata = newdataDF,
                      type = "response",se.fit = TRUE)
etahatschoolg <- predObjschool$fit
lowschoolg <- etahatschoolg - qnorm(0.975)*predObjschool$se.fit
uppschoolg <- etahatschoolg + qnorm(0.975)*predObjschool$se.fit

plot(0,type = "n",xlim = range(OFPforAna$school),ylim = c(2.5,8),
     xlab  =  "number of years of education",
     ylab  =  "mean no. of physician office visits",main  =  "response scale")
polygon(c(schoolg,rev(schoolg)),c(lowschoolg,rev(uppschoolg)),col = "palegreen",border = FALSE)
lines(schoolg,etahatschoolg,col = "darkgreen",lwd = 2)
rug(jitter(OFPforAna$school),col = "dodgerblue",quiet = TRUE)

############ End of OFPGAMfit ############
