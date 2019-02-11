########## R script: SydneyDisplaySophis ##########

# For doing a sophisticated display of the geographical
# component of the geoadditive model fit to the
# Sydney real estate data.

# Last changed: 25 FEB 2018

# Load required packages:

library(fields) ; library(HRW) ; library(mgcv)

# Extract various variables from the Sydney Real Estate
# data frame:  

data(SydneyRealEstate)
logSalePrice <- SydneyRealEstate$logSalePrice 
longitude <- SydneyRealEstate$longitude
latitude <- SydneyRealEstate$latitude
income <- SydneyRealEstate$income
PM10 <- SydneyRealEstate$PM10

# Fit a geoadditive model:

fitGeoadd <- gam(logSalePrice ~ s(income,k = 15) + s(PM10) 
                  + s(longitude,latitude,bs = "tp",k = 80),
                      method = "REML",data=SydneyRealEstate)

# Set up a plotting mesh:

ngrid <- 201
longrid <- seq(150.5,151.35,length = ngrid)
latgrid <- seq(-34.2,-33.52,length = ngrid)
lonlatMesh <- as.matrix(expand.grid(longrid,latgrid))
lonlatMesh <- as.data.frame(lonlatMesh)
lonlatMeshPlus <- data.frame(PM10 = mean(PM10),
                  income = mean(income),lonlatMesh)
names(lonlatMeshPlus) <- c("PM10","income",
                           "longitude","latitude")
fitmesh <- predict(fitGeoadd,lonlatMeshPlus)

# Load the boundary polygon stored in the 
# `SydneyRealEstateBdry' object:

data(SydneyRealEstateBdry)

# Switch off the pixels outside the boundary polygon:

outInds <- (1:ngrid^2)[pointsInPoly(lonlatMesh, 
   SydneyRealEstateBdry) == FALSE] 
fitmesh[outInds] <- NA 
fitmesh <- matrix(fitmesh,ngrid,ngrid)

# Do the image plot with restriction to the polygon:

image.plot(longrid,latgrid,fitmesh,bty = "l",
      xlab = "degrees longitude", 
      ylab = "degrees latitude",col = terrain.colors(1000),
      legend.args = list(text="effect on mean log(sale price)",
      cex = 1.5,adj = 0.8),axis.args = list(cex.axis = 1.5),
      cex.lab = 1.5,cex.axis = 1.5,
      xlim = c(150.6,151.45),ylim = c(-34.15,-33.55))
lines(SydneyRealEstateBdry,lwd = 2)

# Add locations and names of 15 Sydney suburbs:

suburbNames <- c("Blacktown","Campbelltown","Coogee",
                 "Cronulla","Dee Why","Gordon",
                 "Hornsby","Hunter's Hill","Hurstville",
                 "Liverpool","Palm Beach","Parramatta",
                 "Penrith","Strathfield","Vaucluse")
suburbsLonLatMat <- rbind(c(150.9063,-33.7710),
     c(150.8142,-34.0650),c(151.2555,-33.9190),
     c(151.1522,-34.0574),c(151.2854,-33.7544),
     c(151.1492,-33.7573),c(151.0990,-33.7049),
     c(151.1437,-33.8337),c(151.1000,-33.9667),
     c(150.9231,-33.9209),c(151.3217,-33.6011),
     c(151.0011,-33.8150),c(150.7000,-33.7500),
     c(151.0831,-33.8808),c(151.2712,-33.8558))
for (isub in 1:nrow(suburbsLonLatMat))
{
   points(suburbsLonLatMat[isub,1],suburbsLonLatMat[isub,2],
          col = "navy",pch = 16)
   text(suburbsLonLatMat[isub,1] + 0.01,
        suburbsLonLatMat[isub,2],
        suburbNames[isub],adj = 0,cex = 0.9)
}

############ End of SydneyDisplaySophis ############

