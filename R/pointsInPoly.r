########## R-function: pointsInPoly ##########

# For determining those points are inside
# a polygon specified by "polygon".

# Last changed: 22 FEB 2018

pointsInPoly <- function(pointsCoords,polygon)
{
   # Make first and last vertices equal if not
   # already the case:

   if (polygon[1,1]!=polygon[nrow(polygon),1]) 
      polygon <- rbind(polygon,polygon[,1])
   if (polygon[1,2]!=polygon[nrow(polygon),2]) 
      polygon <- rbind(polygon,polygon[,1])

   numVerts <- nrow(polygon) - 1

   # Translate to coords centered at each point:

   x <- -outer(pointsCoords[,1],polygon[1:numVerts,1],"-")
   y <- -outer(pointsCoords[,2],polygon[1:numVerts,2],"-") 

   # Use tan(theta) = (v1 x v2) / (v1 * v2):

   i2 <- c(2:numVerts,1)
   arg1 <- x * y[,i2] - y * x[,i2]
   arg2 <-  x * x[,i2] + y * y[,i2]
   theta <- atan2(arg1,arg2)

   # Sum angles about each point:

   theta <- abs(theta %*% rep(1,numVerts))

   # Return indicators of whether the inputted points
   # are inside or outside the polygon:

   return(as.vector(ifelse(abs(theta - 2*pi)<1e-06,TRUE,FALSE)))
}

########## End of pointsInPoly ##########
