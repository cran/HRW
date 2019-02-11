\name{growthIndiana}
\alias{growthIndiana}
\usage{data(growthIndiana)}
\title{Adolescent somatic growth in Indiana, U.S.A.}
\description{Data on adolescent somatic growth obtained from a study of the mechanisms of human hypertension development conducted at the Indiana University School of Medicine, Indianapolis, Indiana, U.S.A. The data are restricted to a subset of 216 adolescents in the original study who had at least 9 height measurements. There are a total of 4,123 height measurements taken approximately every 6 months.
}
\format{
  A data frame with 4,123 observations on the following 5 variables:
  \describe{
    \item{\code{idnum}}{identification numbers of the 216 adolescents.}
    \item{\code{height}}{height in centimeters.}
    \item{\code{age}}{age in years.}
    \item{\code{male}}{indicator of the adolescent being male:\cr
                        1 = adolescent is male,\cr
                        0 = adolescent is female.}
    \item{\code{black}}{indicator of the adolescent being black:\cr
                        1 = adolescent is black,\cr
                        0 = adolescent is not black.}
  }
}
\references{Pratt, J.H., Jones, J.J., Miller, J.Z., Wagner, M.A. and Fineberg, N.S. (1989). Racial differences in aldosterone excretion and plasma aldosterone concentrations in children. \emph{New England Journal of Medicine}, \bold{321}, 1152-1157.}
\examples{
library(HRW) ; data(growthIndiana)
growthINblackMales <- growthIndiana[(growthIndiana$male == 1) & (growthIndiana$black == 1),]
if (require("lattice"))
   xyplot(height ~ age|factor(idnum),data = growthINblackMales)
}
