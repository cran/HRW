\name{createBoundary}
\alias{createBoundary}
\usage{createBoundary(x,y)}
\title{Boundary polygon creation}
\description{
Create a boundary polygon corresponding nominally to the effective probability density support of a bivariate dataset via an interactive graphical interface and mouse (or, possibly, touchpad) posititionings and button clicks.
}
\arguments{
  \item{x}{vector containing the x-coordinates of a bivariate dataset.}
  \item{y}{vector containing the y-coordinates of a bivariate dataset.}
}
\value{A two-column matrix containing the vertices of the selected
boundary polygon.}
\details{After the bivariate dataset is displayed on the screen a boundary polygon is selected by performing left mouse (or, possibly, touchpad) clicks on the screen to specify vertex positions, and then moving around in a clockwise direction until the polygon is completed. Completion is achieved by clicking inside the red octagon surrounding the starting vertex.}
\author{M.P. Wand 
\email{matt.wand@uts.edu.au}
}
\seealso{
\code{\link{pointsInPoly}}
}
\examples{
library(HRW)
x <- c(4,1,9,8,3,9,7)
y <- c(5,7,5,4,2,1,1)
\dontrun{myBoundary <- createBoundary(x,y)}
}