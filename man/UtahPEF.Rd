\name{UtahPEF}
\alias{UtahPEF}
\usage{data(UtahPEF)}
\title{Peak expiratory flow in Utah, U.S.A.}
\description{
The \code{UtahPEF} data frame data contains longitudinal data on peak expiratory flow, air pollution and temperature for a cohort of 41 schoolchildren in the Utah Valley, U.S.A., during 1990-1991.
}
\format{
This data frame contains the following columns:
\describe{
       \item{\code{idnum}}{schoolchild identification number.}
       \item{\code{devPEF}}{daily peak expiratory flow measurements for each schoolchild minus the overall average for that schoolchild.}
       \item{\code{PM10withMA5}}{5-day moving average of the concentration of particulate matter 10 micrometers or less in diameter.}
       \item{\code{lowTemp}}{lowest temperature in degrees Fahrenheit on the day of recording.}
       \item{\code{timeTrend}}{day of the study on which the measurements were made.}
}
}
\source{Pope, C.A., Dockery, D.W., Spengler, J.D. and Raizenne, M.E. (1991). Respiratory health and PM_10 pollution: a daily time series analysis. \emph{American Review of Respiratory Disease}, \bold{144}, 668-674.
}
\examples{
library(HRW) ; data(UtahPEF)
if (require("lattice"))
   print(xyplot(devPEF ~ PM10withMA5|as.factor(idnum),data = UtahPEF))
}
