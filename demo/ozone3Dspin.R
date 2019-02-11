########## R script: ozone3Dspin ##########

# For illustrating three-dimension spin graphics
# display of a bivariate spline fit.

# Last changed: 07 SEP 2016 

# Load required packages:

library(HRW) ; library(rgl) ; library(mgcv)

# Set function clean():

clean <- function(x) return(x[is.na(x)==0])

# Set up data:

data(ozoneSub)
x1 <- ozoneSub$longitude
x2 <- ozoneSub$latitude
y <- ozoneSub$ozone

# Obtain bivariate spline fit:  

fit <- gam(y ~ s(x1,x2,k = 60))

# Obtain predicting values over a rectangular mesh:

ngrid <- 201
x1grid <- seq(min(x1),max(x1),length = ngrid)
x2grid <- seq(min(x2),max(x2),length = ngrid)
x1x2mesh <- expand.grid(x1grid,x2grid)
names(x1x2mesh) <- c("x1","x2")
fitmesh <- matrix(predict(fit,newdata = x1x2mesh),ngrid,ngrid)

# Set a boundary polygon (obtained using HRW:::createBoundary()):

lonlatBdry <- rbind(c(-87.36,44.58),
                    c(-83.25,43.27),
                    c(-82.79,42.52),
                    c(-82.70,39.89),
                    c(-84.39,37.93),
                    c(-85.99,37.07),
                    c(-88.10,36.64),
                    c(-93.42,37.08),
                    c(-93.71,41.63),
                    c(-91.29,43.93),
                    c(-87.36,44.58))
outInds <- (1:ngrid^2)[pointsInPoly(x1x2mesh,
                       lonlatBdry)==FALSE]
fitmesh[outInds] <- NA

# Plot the fit using three-dimensional spin graphics via
# the package rgl(). For this we work linearly transformed
# versions of the data that are within the unit cube: 

x1UC <- (x1 - min(x1))/(max(x1) - min(x1))
x2UC <- (x2 - min(x2))/(max(x2) - min(x2))
yUC <- (y - min(y))/(max(y) - min(y))

x1gridUC <- (x1grid - min(x1))/(max(x1) - min(x1))
x2gridUC <- (x2grid - min(x2))/(max(x2) - min(x2))
fitmeshUC <- (fitmesh - min(y))/(max(y) - min(y))

# Set up `Right-Hand Coordinate System' versions
# of required rgl functions, which are more natural
# according to the author of this script:

rhcs.spheres <- function(x,y,z,radius,...)
   return(rgl.spheres(x,z,-y,radius,...))

rhcs.texts <- function(x,y,z,text,...)
   return(rgl.texts(x,z,-y,text,...))

rhcs.quads <- function(x,y,z,...)
   return(rgl.quads(x,z,-y,...))

rhcs.surface <- function(x,y,z,...)
   return(rgl.surface(x,z,-y,...))

rhcs.clear <- function(type = "shapes")
{
   rgl.viewpoint(theta = 135)
   return(rgl.clear(type))
}

rhcs.lines <- function (x,y,z,add = FALSE, ...) 
   return(rgl.lines(x,z,-y,...))

# Set colour, alpha-level, point size and axis variables:

bg.col <- "white"     ; bs.col <- "ivory"
ax.col <- "slateblue" ; tx.col <- "navy"      
al.val <- 0.3         ; sp.rad <- 0.015
xtx.pos <- 1.2 ; ytx.pos <- 1.2 ; ztx.pos <- 1.2 
xax.low <- -0.1; yax.low <- -0.1; zax.low <- -0.1
xax.upp <- 1.1 ; yax.upp <- 1.1 ; zax.upp <- 1.1

# Draw the axes and base:   

rhcs.clear()
rgl.bg(col = bg.col)
rhcs.lines(c(xax.low,xax.upp),rep(0,2),rep(0,2),
             size = 3,col = ax.col,add = TRUE)
rhcs.lines(rep(0,2),c(yax.low,yax.upp),rep(0,2),
            size = 3,col = ax.col,add = TRUE)
rhcs.lines(rep(0,2),rep(0,2),c(zax.low,zax.upp),
             size = 3,col = ax.col,add = TRUE)
rhcs.texts(xtx.pos,0,0,"longitude",col = tx.col)
rhcs.texts(0,ytx.pos,0,"latitude",col = tx.col)
rhcs.texts(0,0,ztx.pos,"ozone",col = tx.col)
rhcs.quads(c(xax.low,xax.low,xax.upp,xax.upp),
           c(yax.low,yax.upp,yax.upp,yax.low),
           rep(0,4),col = bs.col,alpha = al.val)
  
# Add the data:

rhcs.spheres(x1UC,x2UC,yUC,col = "dodgerblue",radius = sp.rad)

# Add the surface with shading according to terrain.colors():

al.val <- 0.7
alvec <- rep(al.val,length(fitmeshUC))
numCol <- 501
ycolgrid <- seq(min(clean(fitmeshUC)),max(clean(fitmeshUC)),length = numCol)
ty <- round((numCol-1)*fitmeshUC) +1
colorlut <- c(terrain.colors(numCol))
colvec <- colorlut[ty] 
rgl.surface(x1gridUC,-x2gridUC,fitmeshUC,color = colvec,alpha = alvec)

# Add some points and names of some cities in the region:

Chicagox1 <- -87.6298 ; Chicagox2 <- 41.8781 
Chicagox1UC <- (Chicagox1 - min(x1))/(max(x1) - min(x1))
Chicagox2UC <- (Chicagox2 - min(x2))/(max(x2) - min(x2))
rhcs.spheres(Chicagox1UC,Chicagox2UC,0,col = "darkgreen",radius = 0.5*sp.rad)
rhcs.texts(Chicagox1UC,Chicagox2UC,0.05,"Chicago",col = "darkgreen",cex=0.8)
indx1 <- length(x1gridUC[x1gridUC<=Chicagox1UC])
indx2 <- length(x2gridUC[x2gridUC<=Chicagox2UC])
rhcs.lines(rep(Chicagox1UC,2),rep(Chicagox2UC,2),c(0,fitmeshUC[indx1,indx2]),
             size = 1,col = "darkgreen",add = TRUE)
rhcs.spheres(Chicagox1UC,Chicagox2UC,fitmeshUC[indx1,indx2],col = "darkgreen",radius = 0.5*sp.rad)

Indianapolisx1 <- -86.1581 ; Indianapolisx2 <- 39.7694 
Indianapolisx1UC <- (Indianapolisx1 - min(x1))/(max(x1) - min(x1))
Indianapolisx2UC <- (Indianapolisx2 - min(x2))/(max(x2) - min(x2))
rhcs.spheres(Indianapolisx1UC,Indianapolisx2UC,0,col = "darkgreen",radius = 0.5*sp.rad)
rhcs.texts(Indianapolisx1UC,Indianapolisx2UC,0.05,"Indianapolis",col = "darkgreen",cex=0.8)
indx1 <- length(x1gridUC[x1gridUC<=Indianapolisx1UC])
indx2 <- length(x2gridUC[x2gridUC<=Indianapolisx2UC])
rhcs.lines(rep(Indianapolisx1UC,2),rep(Indianapolisx2UC,2),c(0,fitmeshUC[indx1,indx2]),
             size = 1,col = "darkgreen",add = TRUE)
rhcs.spheres(Indianapolisx1UC,Indianapolisx2UC,fitmeshUC[indx1,indx2],col = "darkgreen",radius = 0.5*sp.rad)

Milwaukeex1 <- -87.9065 ; Milwaukeex2 <- 43.0389
Milwaukeex1UC <- (Milwaukeex1 - min(x1))/(max(x1) - min(x1))
Milwaukeex2UC <- (Milwaukeex2 - min(x2))/(max(x2) - min(x2))
rhcs.spheres(Milwaukeex1UC,Milwaukeex2UC,0,col = "darkgreen",radius = 0.5*sp.rad)
rhcs.texts(Milwaukeex1UC,Milwaukeex2UC,0.05,"Milwaukee",col = "darkgreen",cex=0.8)
indx1 <- length(x1gridUC[x1gridUC<=Milwaukeex1UC])
indx2 <- length(x2gridUC[x2gridUC<=Milwaukeex2UC])
rhcs.lines(rep(Milwaukeex1UC,2),rep(Milwaukeex2UC,2),c(0,fitmeshUC[indx1,indx2]),
             size = 1,col = "darkgreen",add = TRUE)
rhcs.spheres(Milwaukeex1UC,Milwaukeex2UC,fitmeshUC[indx1,indx2],col = "darkgreen",radius = 0.5*sp.rad)

StLouisx1 <- -90.1994 ; StLouisx2 <- 38.6270 
StLouisx1UC <- (StLouisx1 - min(x1))/(max(x1) - min(x1))
StLouisx2UC <- (StLouisx2 - min(x2))/(max(x2) - min(x2))
rhcs.spheres(StLouisx1UC,StLouisx2UC,0,col = "darkgreen",radius = 0.5*sp.rad)
rhcs.texts(StLouisx1UC,StLouisx2UC,0.05,"StLouis",col = "darkgreen",cex=0.8)
indx1 <- length(x1gridUC[x1gridUC<=StLouisx1UC])
indx2 <- length(x2gridUC[x2gridUC<=StLouisx2UC])
rhcs.lines(rep(StLouisx1UC,2),rep(StLouisx2UC,2),c(0,fitmeshUC[indx1,indx2]),
             size = 1,col = "darkgreen",add = TRUE)
rhcs.spheres(StLouisx1UC,StLouisx2UC,fitmeshUC[indx1,indx2],col = "darkgreen",radius = 0.5*sp.rad)

############ End of ozone3Dspin ############
