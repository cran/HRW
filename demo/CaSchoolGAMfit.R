############ R script: CaSchoolGAMfit ##########

# For plotting the estimated function components of
# a generalized additive model fit to the
# California schools data.

# Last changed: 29 NOV 2016

# Load required packages:

library(mgcv) ; library(Ecdat) 

# Load in data:

data(Caschool)
Caschool$log.avginc <- log(Caschool$avginc)

# Obtain GAM fit:

fitGAMCaschool <- gam(mathscr ~ s(calwpct) + s(log.avginc) + 
                                s(compstu) + s(expnstu),data = Caschool)

# Plot the estimated GAM components:

par(mfrow=c(2,2),cex.main = 1.8, col.main = "navy",
    cex.lab = 1.5,mai = c(0.7,0.68,0.1,0.05))
plot(fitGAMCaschool,select = 1,shade = TRUE,bty = "l", 
     xlab = "percent qualifying for CalWORKs",
     ylab = "effect on mean ave. math. score", shade.col = "palegreen",
     rug = FALSE)
rug(Caschool$calwpct,col="dodgerblue")

plot(fitGAMCaschool,select = 2,shade = TRUE,bty = "l", 
     xlab = "log(district average income)",
     ylab = "effect on mean ave. math. score", shade.col = "palegreen",
     rug = FALSE)
rug(Caschool$log.avginc,col="dodgerblue")

plot(fitGAMCaschool,select = 3,shade = TRUE,bty = "l", 
     xlab = "number of computers per student",
     ylab = "effect on mean ave. math. score", shade.col = "palegreen",
     rug = FALSE)
rug(Caschool$compstu,col="dodgerblue")

plot(fitGAMCaschool,select = 4,shade = TRUE,bty = "l", 
     cex.lab = 1.5,xlab = "expenditure per student",
     ylab = "effect on mean ave. math. score", shade.col = "palegreen",
     rug = FALSE)
rug(Caschool$expnstu,col="dodgerblue")

############ End of CaSchoolGAMfit ############


