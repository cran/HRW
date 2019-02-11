\name{BanglaContrac}
\alias{BanglaContrac}
\usage{data(BanglaContrac)}
\title{Contraception use in Bangladesh}
\description{
The \code{BanglaContrac} data frame has multilevel data from the 1988 Bangladesh Fertility Survey. There are data on contraceptive use, number of children, age and urban-dweller status of 1,934 women grouped in 60 districts.
}
\format{
This data frame contains the following columns:
\describe{
       \item{\code{districtID}}{district identification number.}
        \item{\code{usingContraception}}{indicator that woman is using contraception:\cr
                   1 = woman using contraception at time of survey,\cr
                   0 = woman not using contraception at time of survey.}
       \item{\code{childCode}}{numerical code for number of living children:\cr
                   1 = no living children at time of survey,\cr
                   2 = one living child at time of survey,\cr
		   3 = two living children at time of survey,\cr
		   4 = three or more living children at time of survey.}
      \item{\code{ageMinusMean}}{age in years of woman at time of survey, with mean age subtracted}
      \item{\code{isUrban}}{indicator that woman lives in an urban region:\cr
                   1 = woman an urban region dweller at time of survey,\cr
                   0 = woman a rural region dweller at time of survey.}
}
}
\source{
Huq, N.M. and Cleland, J. (1990). \emph{Bangladesh Fertility Survey 1989 (Main Report)}. Dhaka, Bangladesh: National Institute of Population Research and Training.}
\examples{
library(HRW) ; data(BanglaContrac)
if (require("lattice"))
   print(xyplot(jitter(usingContraception) ~ ageMinusMean|factor(districtID),
         groups = childCode, data = BanglaContrac))
}
