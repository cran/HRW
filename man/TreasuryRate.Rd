\name{TreasuryRate}
\alias{TreasuryRate}
\usage{data(TreasuryRate)}
\title{U.S. Treasury rate
}
\description{One-month maturity U.S. Treasury rate during the period 2001-2013.
}
\format{
  A data frame with 3,117 observations on each of the following 2 variables:
  \describe{
    \item{\code{date}}{days from July 31, 2001 until July 10, 2013 with the occasional missing values due to holidays.}
    \item{\code{rate}}{daily one-month maturity U.S. Treasury rate.}
  }
}
\source{
Federal Reserve Bank of St. Louis, U.S.A.
}
\examples{
library(HRW) ; data(TreasuryRate)
TRdate <- as.Date(TreasuryRate$date,"\%m/\%d/\%Y")[!is.na(TreasuryRate$rate)]
TRrate <- TreasuryRate$rate[!is.na(TreasuryRate$rate)]
plot(TRdate,TRrate,type = "l",bty = "l",xlab = "date",ylab = "U.S. Treasury rate")
}

