\name{femSBMD}
\alias{femSBMD}
\usage{data(femSBMD)}
\title{Female spinal bone mineral densities}
\description{
The \code{femSBMD} data frame has longitudinal data on spinal bone mineral, density, age and ethnicity for female youths from a study on bone mineral acquisition. 
}
\format{
  This data frame contains the following columns:
  \describe{
       \item{\code{idnum}}{identification number unique to each subject.}
       \item{\code{spnbmd}}{spinal bone mineral density in grams per square centimeter.}
       \item{\code{age}}{age of subject in years.}
       \item{\code{ethnicity}}{factor corresponding to the subject's
	 ethnicity with levels \code{Asian}, \code{Black}, \code{Hispanic} and \code{White}.}
       \item{\code{black}}{indicator of the subject being black:\cr
                               0 = subject is black\cr
                               1 = subject is not black.}
       \item{\code{hispanic}}{indicator of the subject being Hispanic:\cr
                               0 = subject is Hispanic\cr
                               1 = subject is not Hispanic.}
       \item{\code{white}}{indicator of the subject being white:\cr
                               0 = subject is white\cr
                               1 = subject is not white.}       
     }
     }
\source{
Bachrach, L.K., Hastie, T., Wang, M.-C., Narasimhan, B. and Marcus, R. (1999). Bone mineral acquisition in healthy Asian, Hispanic, Black, and Caucasian youth: a longitudinal study. \emph{Journal of Clinical Endocrinology and Metabolism}, \bold{84}, 4702-4712.
}
\examples{
library(HRW) ; data(femSBMD)
if (require("lattice"))
   print(xyplot(spnbmd ~ age|factor(ethnicity),groups = idnum,
                data = femSBMD,type = "b"))
}

