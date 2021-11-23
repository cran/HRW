############ R script: WarsawAptsGAMfit ##########

# For plotting the estimated function components of
# a generalized additive model fit, with a factor-by-curve 
# interaction, of the Warsaw apartments data.

# Last changed: 02 NOV 2021

# Load required packages:
 
library(mgcv) ; library(HRW) 

# Load in data:

data(WarsawApts)

# Obtain GAM fit:

fit2GAMWarsaw <- gam(areaPerMzloty ~ 
                 s(construction.date,k = 25,by = as.factor(district))
                 + as.factor(n.rooms) + s(surface,k = 25),
                 data = WarsawApts,method = "REML")

# Plot the estimated GAM components:

par(mfrow=c(3,2),mai = c(0.7,0.68,0.3,0.05),cex.lab=1.6,
    cex.axis = 1.6,cex.main = 1.8)

plot(fit2GAMWarsaw,select = 1,shade = TRUE,bty = "l",
     shade.col = "palegreen",col.main = "navy",
     xlab = "construction date (year)",ylab = "",
     main = "Mokotow district",rug = FALSE)
rug(WarsawApts$construction.date,col="dodgerblue")
mtext("effect on mean area",line = 4,side = 2,cex = 1.2)
mtext("per million zloty",line = 3,side = 2,cex = 1.2)

plot(fit2GAMWarsaw,select = 2,shade = TRUE,bty = "l", 
     shade.col = "palegreen", col.main = "navy",
     xlab = "construction date (year)",ylab = "",
     main = "Srodmiescie district",rug = FALSE)
rug(WarsawApts$construction.date,col="dodgerblue")
mtext("effect on mean area",line = 4,side = 2,cex = 1.2)
mtext("per million zloty",line = 3,side = 2,cex = 1.2)

plot(fit2GAMWarsaw,select = 3,shade = TRUE,bty = "l",
     shade.col = "palegreen",col.main = "navy",
     xlab = "construction date (year)",ylab = "",
     main = "Wola district",rug = FALSE)
rug(WarsawApts$construction.date,col="dodgerblue")
mtext("effect on mean area",line = 4,side = 2,cex = 1.2)
mtext("per million zloty",line = 3,side = 2,cex = 1.2)

plot(fit2GAMWarsaw,select = 4,shade = TRUE,bty = "l", 
     shade.col = "palegreen", col.main = "navy",
     xlab = "construction date (year)",ylab = "",
     main = "Zoliborz district",rug = FALSE)
rug(WarsawApts$construction.date,col="dodgerblue")
mtext("effect on mean area",line = 4,side = 2,cex = 1.2)
mtext("per million zloty",line = 3,side = 2,cex = 1.2)

plot(fit2GAMWarsaw,select = 5,shade = TRUE,bty = "l", 
     shade.col = "palegreen", col.main = "navy",
     xlab = "surface area (square meters)",ylab = "",
     rug = FALSE)
rug(WarsawApts$surface,col="dodgerblue")
mtext("effect on mean area",line = 4,side = 2,cex = 1.2)
mtext("per million zloty",line = 3,side = 2,cex = 1.2)

############ End of WarsawAptsGAMfit ############


