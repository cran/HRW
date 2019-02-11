########## R script: MichInctAddMod ##########

# For doing a t-response regression additive
# model example for the "wife work" data
# using the VGAM package.

# Last changed: 19 JUN 2017

# Load required packages:

library(VGAM)  ; library(Ecdat)

# Set-up data:

data(Workinghours)

MichInc <- data.frame(husbandManager = as.numeric(as.character(
                        Workinghours$occupation) == "mp"))
MichInc <- transform(MichInc,
                      otherIncome = Workinghours$income/10,
                      nonWhite = Workinghours$nonwhite,
                      homeOwned = Workinghours$owned,
                      unemployRate = Workinghours$unemp,
                      wifeAge = Workinghours$age,
                      wifeEducationYears = Workinghours$education,
                      numChildren = with(Workinghours, child5 +
                                         child13 + child17))

# Obtain vgam() fit with "family = studentt3":

fit <- vgam(otherIncome ~ s(wifeAge,df = 10)+ s(unemployRate,df = 4)
                        + s(wifeEducationYears,df = 4)
                        + s(numChildren,df = 4) + nonWhite 
                        + homeOwned + husbandManager,
                        family = studentt3,data = MichInc)

# Extract the estimates of the scale and degrees of freedom variables:

sigmaHat <- loge(coef(fit)[2], inverse = TRUE)
nuHat <- loglog(coef(fit)[3], inverse = TRUE)

cat("\n The estimated degrees of freedom for the t-distribution is:",
    signif(nuHat,4),"\n\n")

# Plot the additive fits:

par(mfrow = c(3,3),mai = c(0.7,0.4,0.05,0.04))
plotvgam(fit,se = TRUE,noxmean = TRUE,,bty = "l",lcol = c("darkgreen"),
         scol = "green3",rcol = "dodgerblue",llwd = 2,slwd = 2,ylab = "",
         cex.axis = 1.5,cex.lab = 1.8,scale = 40)

############ End of MichInctAddMod ############


