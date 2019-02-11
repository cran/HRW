\name{pointsInPoly}
\alias{pointsInPoly}
\usage{pointsInPoly(pointsCoords,polygon)}
\title{Points inside/outside polygon determination}
\description{
Determination of whether each member of a set of bivariate points are inside a polygon.
}
\arguments{
  \item{pointsCoords}{two-column matrix with each row specifying the (x,y) coordinates of a point.}
  \item{polygon}{two-column matrix with each row specifying the (x,y) coordinates of a vertex of a polygon.}
}
\value{
A Boolean array with length equal to the number of rows in `pointsCoords' corresponding to whether or not each of the corresponding points in 'pointCoords' are inside 'polygon'.
}
\details{
Geometric results are used to determine whether each of a set of bivariate points is inside or outside a polygon. A Boolean vector of indicators of whether or not each point is inside the polygon is returned.
}
\author{
M.P. Wand \email{matt.wand@uts.edu.au}
}
\seealso{
\code{\link{createBoundary}}
}
\examples{
library(HRW)
myPolygon <- rbind(c(1,9),c(8,8),c(9,3),c(3,2),c(1,9))/10
plot(0:1,0:1,type = "n") ; lines(myPolygon)
xyMat <- cbind(runif(10),runif(10))
inPoly <- pointsInPoly(xyMat,myPolygon) ; print(inPoly)
points(xyMat[,1],xyMat[,2],col = as.numeric(inPoly) + 2)
}
