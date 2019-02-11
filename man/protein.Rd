\name{protein}
\alias{protein}
\usage{data(protein)}
\title{Protein intake dietery study}
\description{
The \code{protein} data frame has longitudinal data on protein intake, body mass index and age of subjects in a dietary study.
}
\format{
  This data frame contains the following columns:
  \describe{
       \item{\code{idnum}}{identification number unique to each subject.}
       \item{\code{proteinBioM}}{logarithm of intake of protein as measured by the biomarker urinary.}
       \item{\code{age}}{age of subject in years.}
       \item{\code{BMI}}{body mass index.}
       \item{\code{proteinRecall}}{logarithm of intake of protein as measured by a 24-hour recall instrument.}
       \item{\code{female}}{indicator that subject is female:\cr
                   1=subject is female,\cr
                   0=subject is male.}}
  }
\source{
Kipnis, V., Subar, A.F., Midthune, D., Freedman, L.S., Ballard-Barbash, R., Troiano, R., Bingham, S., Schoeller, D.A., Schatzkin, A. and Carroll, R.J. (2003). The structure of dietary measurement error: results of the OPEN biomarker study. \emph{American Journal of Epidemiology}, \bold{158}, 14-21.
}
\examples{
library(HRW) ; data(protein)
if (require("lattice"))
   print(xyplot(proteinBioM ~ BMI|factor(female),groups = idnum,
                data = protein, type = "b"))
}
