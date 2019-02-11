\name{coral}
\alias{coral}
\usage{data(coral)}
\title{Coral organisms in French Polynesia}
\description{
The \code{coral} data frame has data on initial size, taxonomic identity and alive/death status of coral organisms in French Polynesia.
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{siteDepthPeriod}}{factor with levels corresponding to a code for the site, depth and time period concerning where and when coral organisms were measured.}
\item{\code{taxon}}{factor corresponding to an abbreviation for taxonomic identity:\cr
  \code{ACR} = Acropora,\cr
  \code{POC} = Pocillopora,\cr
  \code{POR} = Porites.}
\item{\code{logInitialSizePlusOne}}{initial size measurement of coral organism transformed according to the log(initial size + 1).}
\item{\code{died}}{indicator that coral organism has died:\cr
	           1 = coral organism has died,\cr
	           0 = coral organism still alive.
             }
}
}
\source{
Kayal, M., Vercelloni, J., Wand, M.P. and Adjeroud, M. (2015). Searching for the best bet in life-strategy: a quantitative population dynamics approach to life history trade-offs in reef-building corals. \emph{Ecological Complexity}, \bold{23}, 73-84.
}
\examples{
library(HRW) ; data(coral)
if (require("lattice"))
   print(xyplot(died ~ logInitialSizePlusOne|siteDepthPeriod*taxon,
                data = coral,layout = c(15,5)))
}

