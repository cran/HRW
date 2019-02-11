\name{schoolResults}
\alias{schoolResults}
\usage{data(schoolResults)}
\title{School results in the United Kingdom}
\description{The \code{schoolResults} data frame has multilevel data school results and gender for 1,905 school children from 73 schools in United Kingdom.
}
\format{
This data frame contains the following columns:
\describe{
       \item{\code{schoolID}}{school identification number.}
       \item{\code{studentID}}{student identification number.}
       \item{female}{indicator that child is female:\cr
                     1=child is female,\cr
                     0=child is male.}
       \item{\code{writtenScore}}{score on traditional written examination paper out of a total of 160.}
       \item{\code{courseScore}}{score from projects undertaken during the course and marked by the student's own teacher, out of a total of 108.}
}
}  
\source{Creswell, M. (1991). A multilevel bivariate model. In \emph{Data Analysis with ML3} (eds. Prosser, R., Rasbash, J. and Goldstein, H.) London: Institute of Education, London.}
\examples{
library(HRW) ; data(schoolResults)
if (require("lattice"))
   print(xyplot(writtenScore ~ courseScore|factor(schoolID),
                groups = female,data = schoolResults))
}
