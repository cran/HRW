\name{BostonMortgages}
\alias{BostonMortgages}
\usage{data(BostonMortgages)}
\title{Mortgage applications in Boston}
\description{The full dataset is in the data frame `Hdma` within the R package Ecdat. This data frame is the subset of mortgage applications during the years 1998-1999.
}
\format{
  A data frame with 2,029 observations on the following 9 variables:
  \describe{
    \item{\code{dir}}{ratio of the debt payments to the total income.}
    \item{\code{hir}}{ratio of the housing expenses to the total income.}
    \item{\code{lvr}}{ratio of the loan size to the assessed value of property.}
    \item{\code{ccs}}{a credit score ranging from 1 to 6, where a low value indicates low credit risk.}
    \item{\code{mcs}}{a mortgage credit score from 1 to 6, where a low value indicates low credit risk.}
    \item{\code{pbcr}}{did the applicant have a public bad credit record?: a factor with levels \code{no} and \code{yes}.}
    \item{\code{dmi}}{was the applicant denied mortgage insurance?: a factor with levels \code{no} and \code{yes}.}
    \item{\code{self}}{was the applicant self-employed?: a factor with levels \code{no} and \code{yes}.}
    \item{\code{single}}{is the applicant single?: a factor with levels \code{no} and \code{yes}.}
    \item{\code{uria}}{1989 Massachusetts unemployment rate in the applicant's industry.}
    \item{\code{condominium}}{is the unit a condominium?: a factor with levels \code{no} and \code{yes}.}
    \item{\code{black}}{is the applicant black?: a factor with levels \code{no} and \code{yes}.}
    \item{\code{deny}}{was the mortgage denied?: a factor with levels \code{no} and \code{yes}.}
  }
}
\references{
Munnell, A. H., Tootell, G. M. B., Browne, L. E., and McEneaney, J. (1996).
Mortgage lending in Boston: Interpreting HDMA data,
\emph{American Economic Review}, 25-53.
}
\examples{
library(HRW) ; data(BostonMortgages)
BostonMortgages$denyBinary <- as.numeric(BostonMortgages$deny == "yes")
fit <- glm(denyBinary ~ black + dir + lvr + pbcr + self + single + as.factor(ccs),
                        family = binomial,data = BostonMortgages)
summary(fit)
}
