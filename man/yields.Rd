\name{yields}
\alias{yields}
\usage{data(yields)}
\title{U.S., European and Japanese yield curves}
\description{U.S., European and Japanese yield curves. These are functions of maturity.  This dataset has 91 columns. The first column is the date, columns 2 to 31 are European yields at maturities from 1 to 30 years, columns 32 to 61 are Japanese yields at these maturities, and columns 62 to 91 are U.S. yields at the same maturities.
}
\format{
  A data frame with 1565 observations on the following 91 variables:
  \describe{
    \item{\code{date}}{date when the yield is measured.}
    \item{\code{EU01}}{European yield at a maturity of 1 year.}
    \item{\code{EU02}}{European yield at a maturity of 2 years.}
    \item{\code{EU03}}{European yield at a maturity of 3 years.}
    \item{\code{EU04}}{European yield at a maturity of 4 years.}
    \item{\code{EU05}}{European yield at a maturity of 5 years.}
    \item{\code{EU06}}{European yield at a maturity of 6 years.}
    \item{\code{EU07}}{European yield at a maturity of 7 years.}
    \item{\code{EU08}}{European yield at a maturity of 8 years.}
    \item{\code{EU09}}{European yield at a maturity of 9 years.}
    \item{\code{EU10}}{European yield at a maturity of 10 years.}
    \item{\code{EU11}}{European yield at a maturity of 11 years.}
    \item{\code{EU12}}{European yield at a maturity of 12 years.}
    \item{\code{EU13}}{European yield at a maturity of 13 years.}
    \item{\code{EU14}}{European yield at a maturity of 14 years.}
    \item{\code{EU15}}{European yield at a maturity of 15 years.}
    \item{\code{EU16}}{European yield at a maturity of 16 years.}
    \item{\code{EU17}}{European yield at a maturity of 17 years.}
    \item{\code{EU18}}{European yield at a maturity of 18 years.}
    \item{\code{EU19}}{European yield at a maturity of 19 years.}
    \item{\code{EU20}}{European yield at a maturity of 20 years.}
    \item{\code{EU21}}{European yield at a maturity of 21 years.}
    \item{\code{EU22}}{European yield at a maturity of 22 years.}
    \item{\code{EU23}}{European yield at a maturity of 23 years.}
    \item{\code{EU24}}{European yield at a maturity of 24 years.}
    \item{\code{EU25}}{European yield at a maturity of 25 years.}
    \item{\code{EU26}}{European yield at a maturity of 26 years.}
    \item{\code{EU27}}{European yield at a maturity of 27 years.}
    \item{\code{EU28}}{European yield at a maturity of 28 years.}
    \item{\code{EU29}}{European yield at a maturity of 29 years.}
    \item{\code{EU30}}{European yield at a maturity of 30 years.}
    \item{\code{JP01}}{Japanese yield at a maturity of 1 year.}
    \item{\code{JP02}}{Japanese yield at a maturity of 2 years.}
    \item{\code{JP03}}{Japanese yield at a maturity of 3 years.}
    \item{\code{JP04}}{Japanese yield at a maturity of 4 years.}
    \item{\code{JP05}}{Japanese yield at a maturity of 5 years.}
    \item{\code{JP06}}{Japanese yield at a maturity of 6 years.}
    \item{\code{JP07}}{Japanese yield at a maturity of 7 years.}
    \item{\code{JP08}}{Japanese yield at a maturity of 8 years.}
    \item{\code{JP09}}{Japanese yield at a maturity of 9 years.}
    \item{\code{JP10}}{Japanese yield at a maturity of 10 years.}
    \item{\code{JP11}}{Japanese yield at a maturity of 11 years.}
    \item{\code{JP12}}{Japanese yield at a maturity of 12 years.}
    \item{\code{JP13}}{Japanese yield at a maturity of 13 years.}
    \item{\code{JP14}}{Japanese yield at a maturity of 14 years.}
    \item{\code{JP15}}{Japanese yield at a maturity of 15 years.}
    \item{\code{JP16}}{Japanese yield at a maturity of 16 years.}
    \item{\code{JP17}}{Japanese yield at a maturity of 17 years.}
    \item{\code{JP18}}{Japanese yield at a maturity of 18 years.}
    \item{\code{JP19}}{Japanese yield at a maturity of 19 years.}
    \item{\code{JP20}}{Japanese yield at a maturity of 20 years.}
    \item{\code{JP21}}{Japanese yield at a maturity of 21 years.}
    \item{\code{JP22}}{Japanese yield at a maturity of 22 years.}
    \item{\code{JP23}}{Japanese yield at a maturity of 23 years.}
    \item{\code{JP24}}{Japanese yield at a maturity of 24 years.}
    \item{\code{JP25}}{Japanese yield at a maturity of 25 years.}
    \item{\code{JP26}}{Japanese yield at a maturity of 26 years.}
    \item{\code{JP27}}{Japanese yield at a maturity of 27 years.}
    \item{\code{JP28}}{Japanese yield at a maturity of 28 years.}
    \item{\code{JP29}}{Japanese yield at a maturity of 29 years.}
    \item{\code{JP30}}{Japanese yield at a maturity of 30 years.}
    \item{\code{US01}}{U.S. yield at a maturity of 1 year.}
    \item{\code{US02}}{U.S. yield at a maturity of 2 years.}
    \item{\code{US03}}{U.S. yield at a maturity of 3 years.}
    \item{\code{US04}}{U.S. yield at a maturity of 4 years.}
    \item{\code{US05}}{U.S. yield at a maturity of 5 years.}
    \item{\code{US06}}{U.S. yield at a maturity of 6 years.}
    \item{\code{US07}}{U.S. yield at a maturity of 7 years.}
    \item{\code{US08}}{U.S. yield at a maturity of 8 years.}
    \item{\code{US09}}{U.S. yield at a maturity of 9 years.}
    \item{\code{US10}}{U.S. yield at a maturity of 10 years.}
    \item{\code{US11}}{U.S. yield at a maturity of 11 years.}
    \item{\code{US12}}{U.S. yield at a maturity of 12 years.}
    \item{\code{US13}}{U.S. yield at a maturity of 13 years.}
    \item{\code{US14}}{U.S. yield at a maturity of 14 years.}
    \item{\code{US15}}{U.S. yield at a maturity of 15 years.}
    \item{\code{US16}}{U.S. yield at a maturity of 16 years.}
    \item{\code{US17}}{U.S. yield at a maturity of 17 years.}
    \item{\code{US18}}{U.S. yield at a maturity of 18 years.}
    \item{\code{US19}}{U.S. yield at a maturity of 19 years.}
    \item{\code{US20}}{U.S. yield at a maturity of 20 years.}
    \item{\code{US21}}{U.S. yield at a maturity of 21 years.}
    \item{\code{US22}}{U.S. yield at a maturity of 22 years.}
    \item{\code{US23}}{U.S. yield at a maturity of 23 years.}
    \item{\code{US24}}{U.S. yield at a maturity of 24 years.}
    \item{\code{US25}}{U.S. yield at a maturity of 25 years.}
    \item{\code{US26}}{U.S. yield at a maturity of 26 years.}
    \item{\code{US27}}{U.S. yield at a maturity of 27 years.}
    \item{\code{US28}}{U.S. yield at a maturity of 28 years.}
    \item{\code{US29}}{U.S. yield at a maturity of 29 years.}
    \item{\code{US30}}{U.S. yield at a maturity of 30 years.}
  }
}
\source{European Central Bank (European yields), Japanese Ministry of Finance (Japanese yields)
and U.S. Federal Reserve Board (U.S. yields).
}
\examples{
library(HRW) ; data(yields)
t <- 1:30 ; yieldsCleaned <- na.omit(yields)[,-1]
plot(t,yieldsCleaned[1,61:90], type="l",ylim = c(0,6),lwd = 2,
     xlab = "maturity (years)",ylab = "yield",
     bty = "l",cex.lab = 1.5,cex.axis = 1.5)
for (i in 2:14) lines(t,yieldsCleaned[100*i+1,61:90],col = i,lwd = 2)
}

