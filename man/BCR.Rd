\name{BCR}
\alias{BCR}
\usage{data(BCR)}
\title{Berry, Carroll and Ruppert (2002) data}
\description{
The \code{BCR} data frame has data from a 6-week clinical trial of a drug versus a placebo. The data are subject to measurement error and have been transformed and rescaled from their original form. These data are analyzed in the 2002 `Journal of the American Statistical Association' article by Berry, Carroll and Ruppert (full reference below).
}
\format{
This data frame contains the following columns:
\describe{
\item{\code{status}}{code for status of patient:\cr
         control = patient is in placebo group,\cr
	 treatment=patient is in treatment group.}
\item{\code{w}}{physician-assessed score of patient's mental health at baseline.}
\item{\code{y}}{physician-assessed score of patient's mental health at end of the study.}
}
}	   
\references{
Berry, S.M., Carroll, R.J. and Ruppert, D. (2002). Bayesian smoothing and regression splines for measurement error problems.\emph{Journal of the American Statistical Association}, \bold{97}, 160-169.
}
\examples{
library(HRW) ; data(BCR)
if (require("lattice")) 
   print(xyplot(y ~ w|status,data = BCR))
}