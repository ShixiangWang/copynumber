\name{plotGamma}
\alias{plotGamma}
\title{
Plot segmentation results for several values of gamma
}
\description{
Data for one sample on one chromosome is segmented by \code{pcf} for 10 values of gamma, and results are visualized in a multi-grid plot.
}
\usage{
plotGamma(data, pos.unit = "bp", gammaRange = c(10,100), dowins = TRUE, 
          sample = 1, chrom = 1, cv = FALSE, K = 5, cex = 2, col = "grey",
          seg.col="red", ...)
}
\arguments{
  \item{data}{either a data frame or the name of a tab-separated file from which copy number data can be read. The rows of the data frame or file should represent the probes. Column 1 must hold numeric or character chromosome numbers, column 2 the numeric local probe positions, and subsequent column(s) the numeric copy number measurements for one or more samples. The header of copy number columns should give sample IDs.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{gammaRange}{a vector of length two giving the lowest and highest value of gamma to be applied. 10 (approximately) equally spaced values within this range are applied in the pcf-segmentation. Default range is \code{c(10,100)}.}
  \item{dowins}{logical value indicating whether data should be winsorized before running \code{pcf}. Default is TRUE.}
  \item{sample}{an integer indicating which sample is to be segmented. The number should correspond to the sample's place (in order of appearance) in \code{data}. Default is to use the first sample present in the data input.}
  \item{chrom}{a number or character indicating which chromosome is to be segmented. Default is chromosome 1.}
  \item{cv}{logical value indicating whether K-fold cross-validation should be done, see details.}
  \item{K}{the number of folds to use in K-fold cross-validation, default is 5.}
  \item{cex}{size of data points, default is 2.}
  \item{col}{color used to plot data points, default is "grey".}
  \item{seg.col}{color used to plot segments, default is "red".}
  \item{\dots}{other optional parameters to be passed to \code{pcf}.}
}
\details{
Data for one sample and one chromosome is selected, and \code{pcf} is run on this data subset while applying 10 different gamma-values (within the given range). The output is a multi-grid plot with the data shown in the first panel, the segmentation results for the various gammas in the subsequent 10 panels, and the number of segments found for each gamma in the last panel.

If \code{cv = TRUE} a K-fold cross-validation is also performed. For each fold, a random (100/K) per cent of the data are set to be missing, and \code{pcf} is run using the different values of \code{gamma}. The missing probe values are then predicted by the estimated value of their closest non-missing neighbour (see \code{pcf} on this), and the prediction error for this fold is then calculated as the sum of the squared difference between the predicted and the observed values. The process is repeated over the K folds, and the average prediction errors are finally plotted along with the number of segments in the last panel of the plot. The value of gamma for which the minimum prediction error is found is marked by an asterix. Note that such cross-validation tends to favor small values of gamma, and the suitability of the so-called optimal gamma from this procedure should be critically assessed.
}

\value{
  If \code{cv = TRUE} a list containing:
  \item{gamma}{the gamma values applied.}
  \item{pred.error}{the average prediction error for each value of gamma.}
  \item{opt.gamma}{the gamma for which the average prediction error is minimized.}
}


\note{
  This function applies \code{par(fig)}, and is therefore not compatible with other setups for arranging multiple plots in one device such as \code{par(mfrow,mfcol)}. 
}

\seealso{
\code{\link{pcf}},\code{\link{winsorize}}
}
\author{
Gro Nilsen, Knut Liestoel, Ole Christian Lingjaerde
}
\examples{
#Micma data
data(micma)

plotGamma(micma,chrom=17)
}