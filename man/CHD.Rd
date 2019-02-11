\name{CHD}
\alias{CHD}
\usage{data(CHD)}
\title{Coronory heart disease}
\description{The \code{CHD} data frame has data on coronary heart disease status, cholesterol level measurements and age. Further details are given in the 1996 `Journal of the American Statistical Association' article by Roeder, Carroll and Lindsay (full reference below).
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{CHD}}{indicator of coronary heart disease status:\cr
	 0 = patient does not have coronary heart disease,\cr
	 1 = patient has coronary heart disease.}
\item{\code{LDL}}{low density lipoprotein cholesterol level.}
\item{\code{TC}}{total cholesterol level.}
\item{\code{age}}{age of patient in years.}
}
}
\source{Roeder, K., Carroll, R.J. and Lindsay, B.G. (1996). A semiparametric mixture approach to case-control studies with errors in covariables. \emph{Journal of the American Statistical Association}, \bold{91}, 722-732.}	 
\references{
Richardson, S., Leblond, L., Jaussent, I. and Green, P.J. (2002). Mixture models in measurement error problems, with reference to epidemiological studies. \emph{Journal of the Royal Statistical Society, Series A}, \bold{163}, 549-566.
}
\examples{
library(HRW) ; data(CHD)
pairs(CHD)
}
