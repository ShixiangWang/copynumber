\name{multipcf}
\alias{multipcf}
\title{Multi-sample copy number segmentation.}
\description{
  Joint segmentation resulting in piecewise constant curves with common break points for all samples.
}
\usage{
multipcf(data, pos.unit = "bp", arms = NULL, Y = NULL, gamma = 40, 
          normalize=TRUE, w=1, fast = TRUE, assembly = "hg19", digits = 4,
          return.est = FALSE, save.res = FALSE, file.names = NULL, verbose 
          = TRUE)
}
\arguments{
  \item{data}{either a data frame or the name of a tab-separated file from which copy number data can be read. The rows of the data frame or file should represent the probes. Column 1 must hold numeric or character chromosome numbers, column 2 the numeric local probe positions, and subsequent columns the numeric copy number measurements for two or more samples. The header of copy number columns should give sample IDs.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{arms}{optional character vector containing chromosome arms (denoted 'p' and 'q') corresponding to the chromosomes and positions found in \code{data}. If not specified chromosome arms are found using the built-in genome assembly version determined by \code{assembly}.}
  \item{Y}{either a data frame or the name of a tab-separated file containing original copy number data in the case where \code{data} contains Winsorized values. If provided, these values are used to calculate the mean of each segment, otherwise the copy number values in \code{data} are used. \code{Y} must be on the same form as \code{data}.}
  \item{gamma}{penalty for each discontinuity in the curve, default is 40.}
  \item{normalize}{a logical value indicating whether each sample's copy number measurements should be scaled by the sample specific residual standard error. Default is TRUE.}
  \item{w}{a numeric vector giving an individual weight to be used for each sample. May be of length 1 if the same weight should be applied for each sample, default is 1 (no weighting).}
  \item{fast}{a logical value indicating whether a fast (not guaranteed to be exact) version should be run on chromosome arms with > 400 probes. Default is TRUE.}
  \item{assembly}{a string specifying which genome assembly version should be applied to determine chromosome arms. Allowed options are "hg19", "hg18", "hg17" and "hg16" (corresponding to the four latest human genome annotations in the UCSC genome browser).}
  \item{digits}{the number of decimals to be applied when reporting results. Default is 4.}
  \item{return.est}{logical value indicating whether a data frame with copy number estimates (multipcf estimates)should be returned along with the segments. Default is FALSE, which means that only segments are returned.}
  \item{save.res}{logical value indicating whether results should be saved in text files, default is FALSE.}
  \item{file.names}{optional character vector of length two giving the name of the files where the multipcf estimates and segments, respectively, should be saved in case \code{save.res = TRUE}.}
  \item{verbose}{logical value indicating whether or not to print a progress message each time multipcf analysis is finished for a new chromosome arm.}
}
\details{
Piecewise constant curves are simultaneously fitted to the copy number data for several samples as described in the multiPCF algorithm in Nilsen and Liestoel et al. (2012). This implies that break points will be the same for all segmentation curves, but the mean segment values will differ among samples. Segmentation is done separately on each chromosome arm. 

}
\value{
  If \code{return.est = TRUE} a list with the following components:
  \item{estimates}{a data frame where the first two columns give the chromosome numbers and probe positions, respectively, while subsequent columns give the copy number estimates for each sample. The estimate for a given probe and sample equals the sample mean of the segment where the probe is located.}
  \item{segments}{a data frame describing the segments found in the data. Each row represents a segment, and the first five columns give the chromosome numbers,
  arms, local start positions, local end positions, and the number of probes in the segments, respectively. Subsequent columns give the mean segment value for each sample, with sample IDs as column headers.}
  
  If \code{return.est = FALSE} only the data frame containing the segments is returned. 
  
  If \code{save.res = TRUE} the results are also saved in text files with names as specified in \code{file.names}. If \code{file.names=NULL}, a folder named "multipcf_results" is created in the working directory, and the segments and copy number estimates are saved in this folder as tab-separated files named segments.txt and estimates.txt, respectively.
  
}

\note{
It is usually advisable to Winsorize data before running pcf, see \code{\link{winsorize}} on this.

The input data must be complete, see \code{\link{imputeMissing}} for imputation of missing copy number values.
}
\seealso{
\code{\link{imputeMissing}}, \code{\link{pcf}}
}

\references{Nilsen and Liestoel et al., "Copynumber: Efficient algorithms for single- and multi-track copy number segmentation", BMC Genomics 13:591 (2012), doi:10.1186/1471-2164-13-59}

\author{Gro Nilsen, Knut Liestoel}

\examples{
#Load lymphoma data:
data(lymphoma)

#Take out a subset of 3 biopsies from the first patient (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:3)

#Check for missing values in data:
any(is.na(sub.lymphoma))
#FALSE

#First winsorize data to handle outliers:
wins.lymph <- winsorize(sub.lymphoma)

#Run multipcf on subset lymphoma data (using a low gamma because of low-density data)
multi.segments <- multipcf(data=wins.lymph,gamma=12,Y=sub.lymphoma)

}