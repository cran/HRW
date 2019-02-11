########### R function: createBoundary ############

# For formation of polygons around data points.

# Last changed: 21 SEP 2016

createBoundary <- function(x,y)
{ 
   # Check inputs:

   if (missing(x)&missing(y))
      stop("Need to input x and y vectors.")

   # Obtain boundary:

   par(mfrow=c(1,1),xpd = TRUE)
   xlimVal <- c(1.1*min(x) - 0.1*max(x),1.1*max(x) - 0.1*min(x))
   ylimVal <- c(1.1*min(y) - 0.1*max(y),1.1*max(y) - 0.1*min(y))
   plot(x,y,type="n",xaxt="n",yaxt="n",bty="n",xlab="",ylab="",
        xlim = xlimVal,ylim = ylimVal)

   if (length(x)<1000) points(x,y,col = "blue")
   if (length(x)>1000) points(x,y,col = "blue",cex = 0.1)

   cat("\n\n")
   cat("        Using the left mouse button click\n")
   cat("        *clockwise* on the points you want\n")
   cat("        to use to form a polygon.\n\n")
   cat("        The polygon should be closed by\n")
   cat("        clicking inside the red octagon\n")
   cat("        corresponding to the starting vertex.\n\n\n")
   
   finished <- FALSE
   numPoints <- 0 ; bdry <- NULL
   while(!finished)
   {
      locOut <- locator(1,type = "p",pch = 4,col = "green3")
      currPoint <- cbind(locOut$x,locOut$y)
      numPoints <- numPoints + 1

      if (numPoints==1)
      {
         # Draw circle around starting point:

         cex.val <- 0.025
         frame.ratio <- 185/135
         edgx <- cex.val*(max(x)-min(x))
         edgy <- frame.ratio*cex.val*(max(y)-min(y)) 

         circ.x.vert <- currPoint[1] + 0.95*edgx*c(1,1/2,-1/2,-1,-1,-1/2,1/2,1,1)
         circ.y.vert <- currPoint[2] + 0.95*edgy*c(1/2,1,1,1/2,-1/2,-1,-1,-1/2,1/2)
   
         polygon(circ.x.vert,circ.y.vert,density=0,col = "red",lwd = 2)
         points(currPoint,pch=4,col = "green3")
      }

      bdry <- rbind(bdry,currPoint)

      if (numPoints>1)
         lines(c(currPoint[1],bdry[(numPoints-1),1]),
               c(currPoint[2],bdry[(numPoints-1),2]),
               lwd = 2,col = "green3")

      # Test to see if polygon is complete:

      pipOut <- pointsInPoly(currPoint,cbind(circ.x.vert,circ.y.vert))
   
      if ((pipOut==TRUE)&(numPoints>1)) finished <- TRUE
   
   }   

   # Set last vertex of polygon to equal the first

   bdry[nrow(bdry),] <- bdry[1,]

   # Plot data and boundary

   xlimVal <- range(c(xlimVal,bdry[,1]))   
   ylimVal <- range(c(ylimVal,bdry[,2]))

   plot(x,y,type="n",xlim=xlimVal,ylim=ylimVal,xaxt="n",yaxt="n",
        bty="n",xlab="",ylab="")

   if (length(x)<1000) points(x,y,col = "blue")
   if (length(x)>1000) points(x,y,col = "blue",cex = 0.1)

   polygon(bdry,density=0,col = "green3",lwd = 2)

   points(bdry,pch=4,col = "green3")

   # Return selected polygon:
  
   return(bdry)
}

############## End createBoundary ###############
