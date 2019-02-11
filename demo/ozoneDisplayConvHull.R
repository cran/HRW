########## R script: ozoneDisplayConvFull ##########

# For doing a more sophisticated display of the penalized
# thin plate spline fit to the U.S. midwest region ozone 
# data with pixels restricted to the convex hull of the
# geographical data.

# Last changed: 23 FEB 2018

# Load required packages and the U.S. midwest ozone data:

library(mgcv) ; library(HRW) ; data(ozoneSub) ; 

# Obtain thin plate spline fit with 60 basis functions:

fitOzoneThinPlate2 <- gam(ozone ~ s(longitude,latitude,bs = "tp",k = 60), 
                          data = ozoneSub,method = "REML")

# Set up 201 by 201 plotting mesh:

ngrid <- 201
lonGrid <- seq(min(ozoneSub$longitude),
           max(ozoneSub$longitude),length = ngrid)
latGrid <- seq(min(ozoneSub$latitude),max(ozoneSub$latitude),
               length = ngrid)
lonlatMesh <- expand.grid(lonGrid,latGrid)
names(lonlatMesh) <- c("longitude","latitude")

# Obtain the fitted surface over the mesh:

fitMesh <- matrix(predict(fitOzoneThinPlate2,
                  newdata = lonlatMesh),ngrid,ngrid)

# Determine the convex hull of the longitude/latitude
# data and switch off pixels (by setting to NA) that
# are outside the convex hull:

chullInds <- chull(ozoneSub$longitude,ozoneSub$latitude)
chullInds <- c(chullInds,chullInds[1])
ozoneBdry <- cbind(ozoneSub$longitude,
                   ozoneSub$latitude)[chullInds,]
outInds <- (1:ngrid^2)[pointsInPoly(lonlatMesh,
                       ozoneBdry)==FALSE]
fitMesh[outInds] <- NA

# Make the image plot of the convex hull-restricted surface:

par(mai = c(1.02,0.95,0.1,1.2)) 
library(fields)
image.plot(lonGrid,latGrid,fitMesh,col = terrain.colors(1000),
      xlab = "degrees longitude",ylab = "degrees latitude",
      legend.args = list(text = "ozone concentration",
      cex = 1,adj = 0.8),axis.args = list(cex.axis = 1),
      xlim = c(-94.103,-82.429),
      ylim = c(36.408,44.836),bty = "l",
      cex.lab = 1,cex.axis = 1)
lines(ozoneBdry,col = "navy")
points(ozoneSub$longitude,ozoneSub$latitude,
       col = "dodgerblue",cex = 0.5)

# Add U.S.A. state borders for the region being displayed:

US(add = TRUE)

# Add points and names for 4 major cities in the region:

cityNames <- c("Chicago","Indianapolis","Milwaukee","St Louis")
cityLonLatMat <- rbind(c(-87.6298,41.8781),
                       c(-86.1581,39.7694),
                       c(-87.9065,43.0389),
                       c(-90.1994,38.6270))
for (icity in 1:length(cityNames))
{
   points(cityLonLatMat[icity,1],cityLonLatMat[icity,2],
          col = "navy",pch = 16)
   text(cityLonLatMat[icity,1] + 0.15,cityLonLatMat[icity,2],
              cityNames[icity],adj = 0,cex = 1)
}

############ End of ozoneDisplayConvHull ############

