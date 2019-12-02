\name{aspcf}
\alias{aspcf}
\title{Allele-specific copy number segmentation. }
\description{
  Joint segmentation of SNP array data resulting in piecewise constant curves with common break points for copy number data and B-allelle frequency data.
}
\usage{
aspcf(logR, BAF, pos.unit = "bp", arms = NULL, kmin = 5, gamma = 40,
      baf.thres=c(0.1,0.9), skew = 3, assembly= "hg19", digits = 4, 
      return.est = FALSE, save.res = FALSE, file.names=NULL, verbose = TRUE)
}
\arguments{
  \item{logR}{either a data frame or the name of a tab-separated file from which copy number data can be read. The rows of the data frame or file should represent the probes. Column 1 must hold numeric or character chromosome numbers, column 2 the numeric local probe positions, and subsequent columns the numeric copy number measurements for one or more samples. The header of copy number column(s) should give sample ID(s).}
  \item{BAF}{either a data frame or the name of a tab-separated file from which B-allelle frequency data can be read. Must be on the same format and size as \code{logR}, with chromosomes and local probe positions in the two first columns, and numeric BAF-measurements for one or more samples in subsequent columns.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{arms}{optional character vector containing chromosome arms (denoted 'p' and 'q') corresponding to the chromosomes and positions found in \code{logR} and \code{BAF}. If not specified chromosome arms are found using the built-in genome assembly version determined by \code{assembly}.}
  \item{kmin}{minimum number of probes in each segment, default is 5.}
  \item{gamma}{penalty for each discontinuity in the curve, default is 40.}
  \item{baf.thres}{a numeric vector of length two giving the thresholds below and above which BAF probes are considered germline homozygous. Must be in the range 0 to 1, default is 0.1 and 0.9 for the lower and upper limit, respectively.}
  \item{skew}{a numeric value used to determine whether there is allelic skewness (one or two bands) in BAF. Default is 3. The larger the value the further the BAF measurements must be from 0.5 to imply two bands.}
  \item{assembly}{a string specifying which genome assembly version should be applied to determine chromosome arms. Allowed options are "hg19", "hg18", "hg17" and "hg16" (corresponding to the four latest human genome annotations in the UCSC genome browser).}
  \item{digits}{the number of decimals to be applied when reporting results. Default is 4.}
  \item{return.est}{logical value indicating whether a data frame holding LogR estimates should be returned along with the segments. Default is FALSE, which means that only segments are returned.}
  \item{save.res}{logical value indicating whether results should be saved in text files, default is FALSE.}
  \item{file.names}{optional character vector of length two giving the name of the files where the logR estimates and segments, respectively, should be saved in case \code{save.res=TRUE}.}
  \item{verbose}{logical value indicating whether or not to print a progress message each time aspcf analysis is finished for a new chromosome arm.}
}

\details{
Piecewise constant curves are simultaneously fitted to the LogR and BAF data as described in Nilsen and Liestoel et al.(2012). This implies that break points will be the same for the LogR and BAF segmentation curves, while segment values differ.
Segmentation is done separately on each chromosome arm in each sample. 
}

\value{
  If \code{return.est = TRUE} a list with the following components:
  \item{logR_estimates}{a data frame where the first two columns give the chromosome numbers and probe positions, respectively, while subsequent column(s) give the LogR estimates for each sample. The estimate for a given probe equals the mean of the segment where the probe is located.}
  \item{segments}{a data frame describing each segment found. Each row represents a segment, and columns give the sample IDs, chromosome numbers, arms, local start positions, local end positions, number of probes in the segments, mean LogR values and mean BAF values, respectively.}
  
  If \code{return.est = FALSE}, only the data frame containing the segments is returned. 
  
  If \code{save.res = TRUE} the results are also saved in text files with names as specified in \code{file.names}. If \code{file.names=NULL}, a folder named "aspcf_results" is created in the working directory, and the LogR estimates and the segmentation results are saved in this folder as tab-separated files named logR_estimates.txt and segments.txt, respectively.
}

\note{
It will usually be advisable to Winsorize the logR data before running \code{aspcf}, see \code{\link{winsorize}} on this. Missing values are not allowed in \code{logR}, see \code{imputeMissing} for imputation of missing copy number values. 
}

\seealso{\code{\link{plotAllele}}, \code{\link{winsorize}}}
                                                              
\references{Nilsen and Liestoel et al., "Copynumber: Efficient algorithms for single- and multi-track copy number segmentation", BMC Genomics 13:591 (2012), doi:10.1186/1471-2164-13-59}

\author{Gro Nilsen, Knut Liestoel, Ole Christian Lingjaerde}


\examples{
#Load LogR and BAF data:
data(logR)
data(BAF)

#First winsorize logR to handle outliers:
wins.logR <- winsorize(logR)

#Run aspcf:
aspcf.segments <- aspcf(wins.logR,BAF)

}
