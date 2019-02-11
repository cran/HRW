\name{capm}
\alias{capm}
\usage{data(capm)}
\title{Stock indices}
\description{
Daily returns on the Standard & Poors' 500 stock market index, daily rate of the U.S. Treasury bills, and 3 companies' stocks including Microsoft, the General Electric and the Ford Motor Company during the period from November 1, 1993 to March 31, 2003.
}
\format{
  A data frame with 2363 observations on the following 6 variables:
  \describe{
    \item{\code{Close.tbill}}{Daily Treasury bill rate expressed as a percentage.}
    \item{\code{Close.msft}}{Daily closing price of the Microsoft stock.}
    \item{\code{Close.sp500}}{Daily closing Standard and Poor's 500 index.}
    \item{\code{Close.ge}}{Daily closing price of the General Electric stock.}
    \item{\code{Close.ford}}{Daily closing price of the Ford Motor Company stock.}
    \item{\code{Date}}{Dates from November 1, 1993 to March 31, 2003 (d-Mon-yy and dd-Mon-yr formats).}
  }
}
\source{
Federal Reserve Bank of St. Louis U.S.A. (Treasury bill rates) and Yahoo Finance (stock prices).
}
\examples{
# The Capital Asset Pricing Model (CAPM) states that the excess returns on a stock 
# have a linear relationship with the returns on the market. This example investigates
# the CAPM for General Electric stock:

library(HRW) ; data(capm)
n <- dim(capm)[1]
riskfree <- capm$Close.tbill[2:n]/365
elrGE <- diff(log(capm$Close.ge)) - riskfree
elrSP500 <- diff(log(capm$Close.sp500)) - riskfree
plot(elrSP500,elrGE,col = "blue",cex = 0.2)
fitOLS <- lm(elrGE ~ elrSP500)
summary(fitOLS)
par(mfrow = c(2,2)) ; plot(fitOLS)
}

