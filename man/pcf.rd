\name{pcf}
\alias{pcf}
\title{Single-sample copy number segmentation.}
\description{
 Fit a individual piecewise constant segmentation curve to each sample's copy number data.
}
\usage{
pcf(data, pos.unit = "bp", arms = NULL, Y = NULL, kmin = 5, gamma = 40,
      normalize = TRUE, fast = TRUE, assembly = "hg19", digits = 4,
      return.est = FALSE, save.res = FALSE, file.names = NULL, verbose = TRUE)
}

\arguments{
  \item{data}{either a data frame or the name of a tab-separated file from which copy number data can be read. The rows of the data frame or file should represent the probes. Column 1 must hold numeric or character chromosome numbers, column 2 the numeric local probe positions, and subsequent column(s) the numeric copy number measurements for one or more samples. The header of copy number columns should give sample IDs.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{arms}{optional character vector containing chromosome arms (denoted 'p' and 'q') corresponding to the chromosomes and positions found in \code{data}. If not specified chromosome arms are found using the built-in genome assembly version determined by \code{assembly}.}
  \item{Y}{either a data frame or the name of a tab-separated file containing original copy number data in the case where \code{data} contains Winsorized values. If provided, these values are used to calculate the mean of each segment, otherwise the copy number values in \code{data} are used. \code{Y} must be on the same form as \code{data}.}  
  \item{kmin}{minimum number of probes in each segment, default is 5.}
  \item{gamma}{penalty for each discontinuity in the curve, default is 40.}
  \item{normalize}{logical value indicating whether the copy number measurements should be scaled by the sample residual standard error. Default is TRUE.}
  \item{fast}{a logical value indicating whether a fast (not guaranteed to be exact) version should be run on chromosome arms with > 400 probes.}
  \item{assembly}{a string specifying which genome assembly version should be applied to determine chromosome arms. Allowed options are "hg19", "hg18", "hg17" and "hg16" (corresponding to the four latest human genome annotations in the UCSC genome browser).}
  \item{digits}{the number of decimals to be applied when reporting results. Default is 4.}
  \item{return.est}{logical value indicating whether a data frame holding copy number estimates (pcf values) should be returned along with the segments. Default is FALSE, which means that only segments are returned.}
  \item{save.res}{logical value indicating whether results should be saved in text files.}
  \item{file.names}{optional character vector of length two giving the name of the files where the pcf estimates and segments, respectively, should be saved in case \code{save.res=TRUE}.}
  \item{verbose}{logical value indicating whether or not to print a progress message each time pcf analysis is finished for a new chromosome arm.}
}
\details{
  A piecewise constant segmentation curve is fitted to the copy number observations as described in the PCF algorithm in Nilsen and Liestoel et al. (2012). Segmentation is done separately on each chromosome arm in each sample. 
 
}
\value{
  If \code{return.est = TRUE} a list with the following components:
  \item{estimates}{a data frame where the first two columns give the chromosome numbers and probe positions respectively, while subsequent column(s) give the copy number estimates for each sample. The estimate for a given probe equals the mean of the segment where the probe is located.}
  \item{segments}{a data frame describing each segment found in the data. Each row represents a segment, while columns give the sampleID, chromosome number, arm, local start position, local end position, number of probes in the segment and mean value, respectively.}
  
  If \code{return.est = FALSE}, only the data frame containing the segments is returned. 
  
  If \code{save.res = TRUE} the results are also saved in text files with names as specified in \code{file.names}. If \code{file.names=NULL}, a folder named "pcf_results" is created in the working directory, and the pcf estimates and segments are saved in this directory in tab-separated files named estimates.txt and segments.txt, respectively.
}
                                                              

\references{Nilsen and Liestoel et al., "Copynumber: Efficient algorithms for single- and multi-track copy number segmentation", BMC Genomics 13:591 (2012), doi:10.1186/1471-2164-13-59}

\author{Gro Nilsen, Knut Liestoel, Ole Christian Lingjaerde.}

\note{
It is usually advisable to Winsorize data before running pcf, see \code{\link{winsorize}} on this.

Missing copy number values are allowed. These are kept out of the pcf analysis, and copy number estimates for missing observations are later set to be the same as the estimate of the nearest observed probe.
}
\seealso{
\code{\link{multipcf}}
}

\examples{
#Load the lymphoma data set:
data(lymphoma)

#Take out a smaller subset of 3 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:3)

#First winsorize data to handle outliers:
wins.lymph <- winsorize(sub.lymphoma)

#Run pcf (using small gamma because of low-density data):
pcf.segments <- pcf(data=wins.lymph,gamma=12,Y=sub.lymphoma)


}
