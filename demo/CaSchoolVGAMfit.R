############ R script: CaSchoolVGAMfit ##########

# For plotting the estimated function components of a
# vector generalized additive model fit of the California 
# school data with vector response:
#
#     mathscr - average mathematics scores
#     readscr - average reading scores
#
# and five predictors.

# Last changed: 02 JAN 2019

# Load required packages:

library(VGAM) ; library(Ecdat) 

# Load in data:

data(Caschool)
Caschool$log.avginc <- log(Caschool$avginc)

# Obtain VGAM fit:

fitVGAMCaschool <- vgam(cbind(mathscr,readscr) ~ mealpct
   + elpct + s(calwpct,df = 3) + s(compstu,df = 3) 
   + s(log.avginc,df = 4),family = uninormal,data = Caschool)

# Plot the estimated VGAM components:

par(mfrow = c(2,3),lwd = 2,bty = "l",cex.lab = 1.4,
    cex.main = 1.8,cex.axis = 1.4,mai = c(0.7,0.6,0.35,0.05))
curveCols <- c("darkgreen","indianred3")

plot(fitVGAMCaschool,which.term = 1,
     se = TRUE,overlay = TRUE, bty = "l",
     lcol = curveCols, scol = 1:2,  
     xlab="perc. qualif. reduced priced lunch",
     ylab="effect on mean response",rug = FALSE,ylim = c(-25,20))
rug(Caschool$mealpct,col = "dodgerblue")
legend("topright",legend=c("math. score","read. score"),
       col=c("darkgreen","indianred3"),lty=rep(1,2),lwd=rep(2,2))

plot(fitVGAMCaschool,which.term = 2,
     se = TRUE,overlay = TRUE, bty = "l",lcol = curveCols, scol = 1:2, 
     col.main = "navy",xlab="percent of English learners",
     ylab="effect on mean response",rug = FALSE,ylim = c(-25,20))
rug(Caschool$elpct,col = "dodgerblue")

plot(fitVGAMCaschool,which.term = 3,
     se = TRUE,overlay = TRUE, bty = "l",lcol = curveCols, scol = 1:2,
     col.main = "navy",xlab="percent qualifying for CalWORKs",
     ylab="effect on mean response",rug = FALSE,ylim = c(-25,20))
rug(Caschool$calwpct,col = "dodgerblue")

plot(fitVGAMCaschool,which.term = 4,
     se = TRUE,overlay = TRUE, bty = "l",lcol = curveCols, scol = 1:2,
     col.main = "navy",xlab="number of computers per student",
     ylab="effect on mean response",rug = FALSE,ylim = c(-25,20))
rug(Caschool$compstu,col = "dodgerblue")

plot(fitVGAMCaschool,which.term = 5,
     se = TRUE,overlay = TRUE, bty = "l",lcol = curveCols, scol = 1:2,
     col.main = "navy",xlab="log(district average income)",
     ylab="effect on mean response",rug = FALSE,ylim = c(-25,20))
rug(Caschool$log.avginc,col = "dodgerblue")

############ End of CaSchoolVGAMfit ############


