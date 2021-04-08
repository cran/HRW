########## R script: BostMortFacByCurv ##########

# For creation of the Chapter 3 figure of the Harezlak, Ruppert
# and Wand book concerning a factor-by-curve interaction
# model for the Boston Mortgage data.

# Last changed: 06 APR 2021

# Load data:

library(HRW) ; data(BostonMortgages) ; library(mgcv) 
BostonMortgages$denyBinary <- as.numeric(BostonMortgages$deny == "yes")

# Fit model:

fitInt <- gam(denyBinary ~ black + s(dir,by = factor(self),k = 27)
                                 + s(lvr,by = factor(self),k = 27)
                                 + pbcr + self + single + s(ccs,k = 4), 
                                 family = binomial, data = BostonMortgages)

# Set up grids for customised plots of smooth fits:

ng <- 1001
dirLow <- 0 ; dirUpp <- 1
dirg <- seq(dirLow,dirUpp,length = ng)
lvrLow <- 0 ; lvrUpp <- 1
lvrg <- seq(lvrLow,lvrUpp,length = ng)

dirAveg <- rep(mean(BostonMortgages$dir),ng)
lvrAveg <- rep(mean(BostonMortgages$lvr),ng)

selfNog <- as.factor(rep("no",ng))
selfYesg <- as.factor(rep("yes",ng))

blackYesg <- as.factor(rep("yes",ng))
singleYesg <- as.factor(rep("yes",ng))
pbcrYesg <- as.factor(rep("yes",ng))
ccsAveg <- rep(mean(BostonMortgages$ccs),ng)

nonSelfCol <- "indianred3"
selfCol <- "darkgreen"

# Set graphical parameters:

par(mfrow = c(1,2),mai = c(1.02,0.9,0.3,0.2)) ; cex.labVal <- 1.4

# Do plot as a function of "dir":

fdirNog <- predict(fitInt,type = "response",
                   newdata = data.frame(self = selfNog,dir = dirg,
                   lvr = lvrAveg,pbcr = pbcrYesg,black = blackYesg,
                   single = singleYesg,ccs = ccsAveg),se = TRUE)

fdirYesg <- predict(fitInt,type = "response",
                    newdata = data.frame(self = selfYesg,dir = dirg,
                    lvr = lvrAveg,pbcr = pbcrYesg,black = blackYesg,
                    single = singleYesg,ccs = ccsAveg),se = TRUE)

probdirNog <- fdirNog$fit
sdprobNog <- fdirNog$se.fit

lowprobdirNog <- probdirNog - 2*sdprobNog 
uppprobdirNog <- probdirNog + 2*sdprobNog

lowprobdirNog[lowprobdirNog<0] <- 0
uppprobdirNog[uppprobdirNog>1] <- 1

probdirYesg <- fdirYesg$fit
sdprobYesg <- fdirYesg$se.fit

lowprobdirYesg <- probdirYesg - 2*sdprobYesg 
uppprobdirYesg <- probdirYesg + 2*sdprobYesg

lowprobdirYesg[lowprobdirYesg<0] <- 0
uppprobdirYesg[uppprobdirYesg>1] <- 1

plot(0,type = "n",bty = "l",xlim = range(dirg),ylim = c(0,1),
     xlab = "debt to income ratio",ylab = "probability of mortgage application denied",
     cex.lab = cex.labVal)
rug(BostonMortgages$dir,col = "dodgerblue",quiet = TRUE)

lines(dirg,probdirNog,col = nonSelfCol,lwd = 2)
lines(dirg,lowprobdirNog,col = nonSelfCol,lwd = 2,lty = 2)
lines(dirg,uppprobdirNog,col = nonSelfCol,lwd = 2,lty = 2)

lines(dirg,probdirYesg,col = selfCol,lwd = 2)
lines(dirg,lowprobdirYesg,col = selfCol,lwd = 2,lty = 2)
lines(dirg,uppprobdirYesg,col = selfCol,lwd = 2,lty = 2)

abline(h = 0,col = "slateblue",lty = 2)
abline(h = 1,col = "slateblue",lty = 2)

legend(0.42,0.25,legend = c("self employed","not self-employed"),
       lty = rep(1,2),lwd = rep(2,2),col = c(selfCol,nonSelfCol),cex = 0.8)

# Do plot as a function of "dvr":

flvrNog <- predict(fitInt,type = "response",
                   newdata = data.frame(self = selfNog,dir = dirAveg,
                   lvr = lvrg,pbcr = pbcrYesg,black = blackYesg,
                   single = singleYesg,ccs = ccsAveg),se = TRUE)

flvrYesg <- predict(fitInt,type = "response",
                    newdata = data.frame(self = selfYesg,dir = dirAveg,
                    lvr = lvrg,pbcr = pbcrYesg,black = blackYesg,
                    single = singleYesg,ccs = ccsAveg),se = TRUE)

problvrNog <- flvrNog$fit
sdprobNog <- flvrNog$se.fit

lowproblvrNog <- problvrNog - 2*sdprobNog 
uppproblvrNog <- problvrNog + 2*sdprobNog

lowproblvrNog[lowproblvrNog<0] <- 0
uppproblvrNog[uppproblvrNog>1] <- 1

problvrYesg <- flvrYesg$fit
sdprobYesg <- flvrYesg$se.fit

lowproblvrYesg <- problvrYesg - 2*sdprobYesg 
uppproblvrYesg <- problvrYesg + 2*sdprobYesg

lowproblvrYesg[lowproblvrYesg<0] <- 0
uppproblvrYesg[uppproblvrYesg>1] <- 1

plot(0,type = "n",bty = "l",xlim = range(lvrg),ylim = c(0,1),
     xlab = "loan size to property value ratio",
     ylab = "probability of mortgage application denied",
     cex.lab = cex.labVal)
rug(BostonMortgages$lvr,col = "dodgerblue",quiet = TRUE)

lines(lvrg,problvrNog,col = nonSelfCol,lwd = 2)
lines(lvrg,lowproblvrNog,col = nonSelfCol,lwd = 2,lty = 2)
lines(lvrg,uppproblvrNog,col = nonSelfCol,lwd = 2,lty = 2)

lines(lvrg,problvrYesg,col = selfCol,lwd = 2)
lines(lvrg,lowproblvrYesg,col = selfCol,lwd = 2,lty = 2)
lines(lvrg,uppproblvrYesg,col = selfCol,lwd = 2,lty = 2)

abline(h = 0,col = "slateblue",lty = 2)
abline(h = 1,col = "slateblue",lty = 2)

############ End of BostMortFacByCurv ##########
