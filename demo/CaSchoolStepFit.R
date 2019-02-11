############ R script: CaSchoolStepFit ##########

# For plotting the estimated function components of
# a generalized additive model fit of the California schools data,
# with the model selected using a stepwise procedure via the function
# gam:::step.gam().

# Last changed: 08 JAN 2018

# Load required packages:

library(mgcv) ; library(Ecdat)

# Load in data:

data(Caschool)
Caschool$log.avginc <- log(Caschool$avginc)

# Obtain GAM fit:

fitStepCaschool <- gam(mathscr ~   mealpct + elpct
                + s(calwpct,k = 27) + s(compstu,k = 27)
                + s(log.avginc,k =27),data = Caschool)

# Plot the estimated GAM components:

par(mfrow = c(2,3),lwd = 2,bty = "l",cex.lab = 1.4,
    cex.main = 1.8,cex.axis = 1.4,mai = c(0.7,1.2,0.35,0.05))

termplot(fitStepCaschool,terms = c("mealpct","elpct"),se = TRUE,
         col.term = "darkgreen",col.se = "palegreen",rug = TRUE,
         xlabs = c("perc. qualif. reduced priced lunch",
                   "percent of English learners"),
         ylab = rep("effect on mean ave. math. score",2))

plot(fitStepCaschool,select = 1,shade = TRUE,bty = "l",
     xlab = "percent qualifying for CalWORKs",
     ylab = "effect on mean ave. math. score",shade.col = "palegreen",
     rug = FALSE)
rug(Caschool$calwpct,col = "dodgerblue")

plot(fitStepCaschool,select = 2,shade = TRUE,bty = "l",
     xlab = "number of computers per student",
     ylab = "effect on mean ave. math. score",shade.col = "palegreen",
     rug = FALSE)
rug(Caschool$compstu,col = "dodgerblue")

plot(fitStepCaschool,select = 3,shade = TRUE,bty = "l",
     xlab = "log(district average income)",
     ylab = "effect on mean ave. math. score",shade.col = "palegreen",
     rug = FALSE)
rug(Caschool$log.avginc,col = "dodgerblue")

############ End of CaSchoolStepFit ############


