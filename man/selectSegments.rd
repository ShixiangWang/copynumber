\name{selectSegments}
\alias{selectSegments}
\title{
Select multipcf segments
}
\description{
Selects multipcf segments based on a desired characteristic.
}
\usage{
 selectSegments(segments, what = "variance", thres = NULL, nseg = 10, 
                  large = TRUE, p = 0.1)
}
\arguments{
  \item{segments}{a data frame containing segments found by \code{\link{multipcf}}.}
  \item{what}{the desired characteristic to base selection on. Must be one of "variance" (default),"length" and "aberration". See details below.}
  \item{thres}{an optional numeric threshold to be applied in the selection.}
  \item{nseg}{the desired number of segments to be selected, default is 10. Only used if \code{thres=NULL}.}
  \item{large}{logical value indicating whether segments with large (TRUE) or small (FALSE) variance, length or mean value should be selected when \code{what} is "variance", "length" or "aberration", respectively.}
   \item{p}{a number between 0 and 1 giving the minimum proportion of samples for which an aberration must be detected, default is 0.1. Only applicable if \code{what="aberration"}.}
}
\details{
  The input in \code{what} determines how the segments are selected. Three options are available:
  
  If \code{what="variance"} the variance of the segment values across all samples is calculated for each segment. If \code{thres} is specified, the subset of segments for which the variance is above (if \code{large=TRUE}) or below (if \code{large=FALSE}) the threshold is returned. If \code{thres} is not given by the user, a given number of segments determined by the input in \code{nseg} is selected; if \code{large=TRUE} this will be the \code{nseg} segments with the highest variance, whereas if \code{large=FALSE} the subset will consist of the \code{nseg} segments with the lowest variance.
  
  If \code{what="length"} selection is based on the genomic length of the segments (end position minus start position). If \code{thres} is specified, the subset of segments for which the length is above (if \code{large=TRUE}) or below (if \code{large=FALSE}) this threshold is returned. If \code{thres} is left unspecified, a given number of segments determined by the input in \code{nseg} is selected; if \code{large=TRUE} this will be the \code{nseg} longest segments, whereas if \code{large=FALSE} it will be the \code{nseg} shortest segments.
  
  If \code{what="aberration"} the aberration frequency is used to select the subset of segments. If \code{thres} is specified, the proportion of samples for which the segment value is above (if \code{large=TRUE}) or below (if \code{large=FALSE}) the threshold is calculated for each segment. The subset of segments where this frequency is above or equal to the proportion set by the parameter \code{p} is returned. If \code{thres} is not specified, the \code{nseg} segments with the highest (1-p)-quantile (if \code{large=TRUE}) or the lowest p-quantile (if \code{large=FALSE}) is returned.
}
\value{ 
A list containing: 
\item{sel.seg}{data frame containing the selected segments.}
In addition, depending on the value of \code{what}:
\item{seg.var}{a vector giving the variance for each segment. Only returned if \code{what = "variance"}.}
\item{seg.length}{a vector giving the length of each segment. Only returned if \code{what = "length"}.}
\item{seg.ab.prop}{a vector giving the aberration proportion for each segment. Only returned if \code{what = "aberration"} and \code{thres} is specified.}
\item{seg.quantile}{a vector giving the (1-p)- or p-quantile for each segment. Only returned if \code{what = "aberration"} and \code{thres=NULL}.}
}

\seealso{
  \code{\link{multipcf}}
}
\author{Gro Nilsen}

\examples{
#Lymphoma data
data(lymphoma)

#Run multipcf
segments <- multipcf(lymphoma,gamma=12)

#Select the 10 segments with the highest variance:
sel.seg1 <- selectSegments(segments)

#Select the segments where the variance is below 0.001
sel.seg2 <- selectSegments(segments,thres=0.001,large=FALSE)

#Select the 5 longest segments:
sel.seg3 <- selectSegments(segments,what="length",nseg=5)

#Select the segments where 20 % of the samples have segment value of 0.2 or more:
sel.seg4 <- selectSegments(segments,what="aberration",thres=0.2,p=0.2)

#Select the 20 segments with the largest median:
sel.seg5 <- selectSegments(segments,what="aberration",nseg=20,p=0.5)
}
