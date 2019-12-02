\name{callAberrations}
\alias{callAberrations}

\title{
Call aberrations in segmented data
}
\description{
Segments, obtained by \code{pcf} or \code{multipcf}, are classified as "gain", "normal" or "loss" given the specified thresholds.  
}
\usage{
callAberrations(segments, thres.gain, thres.loss = -thres.gain)
}
\arguments{
   \item{segments}{a data frame containing the segmentation results found by either \code{\link{pcf}} or \code{\link{multipcf}}.}
  \item{thres.gain}{a numeric value giving the threshold to be applied for calling gains.}
  \item{thres.loss}{a numeric value giving the threshold to be applied for calling losses. Default is to use the negative value of \code{thres.gain}.}
}
\details{
Each region found in \code{segments} is classified as "gain", "normal" or "loss". Regions with gain or loss will be those segments where the segment value is above or below the value given in \code{thres.gain} or \code{thres.loss}, respectively.
}

\value{
A new segment data frame where the segment values have been replaced by the classification "gain", "normal" or "loss".
}

\author{
Gro Nilsen
}

\examples{
#load lymphoma data
data(lymphoma)
#Run pcf
seg <- pcf(data=lymphoma,gamma=12)

#Call gains as segments whose value is > 0.2, and losses as segments whose
# value < -0.1
ab.seg <- callAberrations(seg,thres.gain=0.2,thres.loss=-0.1)

}

