########## R script: WarsawAptsSimpSemi ##########

# For performing a simple semiparametric regression analysis of
# the Warsaw apartment data with predictor equal to
# year of construction, response equal to the area/price
# ratio and factor equal to district.

# Last changed: 09 JAN 2018

# Load required packages:

library(HRW)  ;  library(mgcv)

# Load Warsaw apartment data:

data(WarsawApts)

# Set graphics parameters:

cex.axisVal <- 1.8 ; cex.labVal <- 1.8
cex.legendVal <- 1.3 ; lwdVal <- 2

# Perform additive factor effects fit and summarise results:

fitSimpSemi <- gam(areaPerMzloty ~ factor(district) +
                     s(construction.date,bs = "cr",k = 27),
                     data = WarsawApts)

print(summary(fitSimpSemi))

# Obtain colour-coded scatterplot with fits:

construcDate <- WarsawApts$construction.date
areaPerMzloty <- WarsawApts$areaPerMzloty
districtChar <- as.character(WarsawApts$district)

par(mai = c(1.02,0.9,0.82,0.42))
plot(construcDate,areaPerMzloty,type="n",bty="l",
     xlab="construction date (year)",
     ylab="area (square meters) per million zloty",
     cex.axis = cex.axisVal,cex.lab = cex.labVal)

myPtCols <- c("deepskyblue","salmon","green3","gold")
myLnCols <- c("blue","red","darkgreen","darkorange")

points(construcDate[districtChar == "Mokotow"],
       areaPerMzloty[districtChar == "Mokotow"],
       lwd = lwdVal,col = myPtCols[1])
points(construcDate[districtChar == "Srodmiescie"],
       areaPerMzloty[districtChar == "Srodmiescie"],
        lwd = lwdVal,col = myPtCols[2])
points(construcDate[districtChar == "Wola"],
       areaPerMzloty[districtChar == "Wola"],
       lwd = lwdVal,col = myPtCols[3])
points(construcDate[districtChar == "Zoliborz"],
       areaPerMzloty[districtChar == "Zoliborz"],
       lwd = lwdVal,col = myPtCols[4])

ng <- 1001
construcDateg <- seq(min(construcDate),max(construcDate),length = ng)
fHatgMokotow <- predict(fitSimpSemi, newdata = data.frame(
                        construction.date = construcDateg,
                        district = rep("Mokotow",ng)))
lines(construcDateg,fHatgMokotow,col = myLnCols[1])
fHatgSrodmiescie <- predict(fitSimpSemi, newdata = data.frame(
                        construction.date = construcDateg,
                        district = rep("Srodmiescie",ng)))
lines(construcDateg,fHatgSrodmiescie,col = myLnCols[2])
fHatgWola <- predict(fitSimpSemi, newdata = data.frame(
                        construction.date = construcDateg,
                        district = rep("Wola",ng)))
lines(construcDateg,fHatgWola,col = myLnCols[3])
fHatgZoliborz <- predict(fitSimpSemi, newdata = data.frame(
                        construction.date = construcDateg,
                        district = rep("Zoliborz",ng)))
lines(construcDateg,fHatgZoliborz,col = myLnCols[4])

legend(1971,80,legend = c("Mokotow","Srodmiescie","Wola","Zoliborz"),
       col = myPtCols,pch = rep(1,4),cex = cex.legendVal,
       pt.cex = rep(1,4),pt.lwd = rep(lwdVal,4))

# Next do a lattice plot with a separate panel for each district:

library(lattice)
tmp <- trellis.par.get("add.text")
tmp$cex <- 1.5
trellis.par.set("add.text",tmp)

pobj <- xyplot(areaPerMzloty ~ construction.date|district,
               data = WarsawApts,as.table=TRUE,
               par.settings = list(layout.heights
                                 =list(strip=1.4)),
               scales = list(cex = 1.25),
               xlab= list("construction date (year)",cex=cex.labVal),
               ylab= list("area (square meters) per million zloty",cex=cex.labVal),
               panel=function(x,y) 
               {
                  panel.grid()
                  panel.xyplot(x,y)
                  if (panel.number() == 1) 
                     panel.xyplot(construcDateg,fHatgMokotow,
                                  col = "darkgreen",lwd = 2,type = "l")
                  if (panel.number() == 2) 
                     panel.xyplot(construcDateg,fHatgSrodmiescie,
                                  col = "darkgreen",lwd = 2,type = "l")
                  if (panel.number() == 3) 
                     panel.xyplot(construcDateg,fHatgWola,
                                  col = "darkgreen",lwd = 2,type = "l")
                  if (panel.number() == 4) 
                     panel.xyplot(construcDateg,fHatgZoliborz,
                                  col = "darkgreen",lwd = 2,type = "l")
               }   
               )
print(pobj)

########## End of WarsawAptsSimpSemi ##########

