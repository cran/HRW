\name{ZOSull}
\alias{ZOSull}
\usage{ZOSull(x,range.x,intKnots,drv = 0)}
\title{O'Sullivan spline design matrices}
\description{Constructs a design matrix consisting of O'Sullivan cubic spline functions of an array of abscissae. Typicially the array corresponds to either observed values of a predictor or an abscissa grid for plotting purposes. 
}
\arguments{
\item{x}{array of abscissae.}
\item{range.x}{array of length 2 such that range.x[1] <= min(x) and range.x[2] >= max(x).}
\item{intKnots}{ordered array of length smaller than that of x and consisting of unique numbers between min(x) and max(x) that specifies the positions of internal knots, that define the spline basis (see the Wand and Ormerod (2008) reference below for full mathematical details).}
\item{drv}{either 0,1 or 2 with a default value of 0. If drv = 1 then the first derivatives of the O'Sullivan spline basis functions are computed instead. Similarly, if drv = 2 then the second derivatives are computed.}
}
\value{A matrix of with length(x) rows and (length(intKnots) + 2) columns with each column containing a separate O'Sullivan spline basis function (determined by range.x and intKnots) of the abscissae in x. The values of range.x and intKnots are included as attributes of the fit object.}
\author{
Matt Wand \email{matt.wand@uts.edu.au}
}
\references{
O'Sullivan, F. (1986). A statistical perspective on ill-posed inverse problems (with discussion). \emph{Statistical Science}, \bold{1}, 505-527.

Wand, M.P. and Ormerod, J.T. (2008). On semiparametric regression with O'Sullivan penalized splines. \emph{Australian and New Zealand Journal of Statistics.} \bold{50}, 179-198. 
}
\examples{
library(HRW)
x <- WarsawApts$construction.date
a <- 1.01*min(x) - 0.01*max(x)
b <- 1.01*max(x) - 0.01*min(x) ; numIntKnots <- 23
intKnots <- quantile(unique(x),seq(0,1,length = (numIntKnots + 2))[-c(1,(numIntKnots + 2))])
xg <- seq(a,b,length = 1001)
Zg <- ZOSull(xg,range.x = c(a,b),intKnots = intKnots)
plot(0,type = "n",xlim = range(xg),ylim = range(Zg),
     bty = "l",xlab = "construction date (year)",
     ylab = "spline basis function")
for (k in 1:ncol(Zg)) lines(xg,Zg[,k],col = k,lwd = 2)
   lines(c(min(xg),max(xg)),rep(0,2),col = "darkmagenta")
for (k in 1:numIntKnots)
   points(intKnots[k],0,pch = 18,cex = 2,col = "darkmagenta")
}




