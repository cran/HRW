\name{brainImage}
\alias{brainImage}
\usage{data(brainImage)}
\title{Brain image}
\description{
The \code{brainImage} data frame corresponds a functional magnetic image of a coronal slice of a human brain.
The data are brain activity on a pixel array.
}
\format{
This data frame is a 80 by 37 array. The entries correspond to brain activity in each of the corresponding pixels. The columns names are c1-c37. These names have no meaning, and are present to ensure that this data frame conforms with R data frame conventions. 
}
\source{
Landman, B.A., Huang, A.J., Gifford, A., Vikram, D.S., Lim, I.A.L, Farrell, J.A.D., Bogovic, J.A., Hua, J., Chen, M., Jarso, S., Smith,  S.A., Joel, S., Mori, S., Pekar, J.J., Barker, P.B., Prince, J.L. and van Zijl, P.C.M. (2010). Multi-parametric neuroimaging reproducibility: A 3T resource study. \emph{NeuroImage}, \bold{54}, 2854-2866. 
}
\examples{
library(HRW) ; data(brainImage) 
image(as.matrix(brainImage))
}
