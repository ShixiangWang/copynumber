\name{interpolate.pcf}
\alias{interpolate.pcf}
\title{Interpolation of pcf-estimates.}
\description{
 Given a segmentation by \code{pcf}, interpolate pcf-estimates for specific positions. 
}
\usage{
interpolate.pcf(segments,x)
}

\arguments{
  \item{segments}{a data frame containing the segmentation results from \code{\link{pcf}}.}
  \item{x}{matrix or data.frame where the first column gives chrosomomes and the second gives positions.}
}
\details{
  Pcf-estimates are interpolated for the chromosomes and postions specified in \code{x}. 
}
\value{
 A data frame where the first two columns give the chromsomes and positions specified in the input \code{x} and the remaining columns give the interpolated pcf-estimate for each sample represented in \code{segments}.
}
                                                              
\author{Gro Nilsen, Ole Christian Lingjaerde.}

\seealso{
\code{\link{pcf}}
}

\note{
The positions in \code{segments} and \code{x} must be of the same unit (bp, kbp, or mbp).
}

\examples{
#Load the lymphoma data set:
data(lymphoma)

#Take out a smaller subset of 3 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:3)

#Run pcf:
seg <- pcf(data=sub.lymphoma,gamma=12)

#Make a matrix with two positions and chromosomes for which we want to 
#interpolate the pcf-estimate:
pos <-  c(2000000,50000000)
chr <- c(1,2)
x <- cbind(chr,pos)

#Interpolate
int.pcf <- interpolate.pcf(seg,x)

}
