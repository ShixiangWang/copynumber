\name{plotGenome}
\alias{plotGenome}

\title{Plot copy number data and/or segmentation results}

\description{                                            
 Plot copy number data and/or segmentation results for the whole genome.
 
}
\usage{
plotGenome(data = NULL, segments = NULL, pos.unit = "bp", sample = NULL,
            assembly="hg19", winsoutliers = NULL, xaxis = "pos", 
            layout = c(1,1), ...)
}
\arguments{
  \item{data}{a data frame with numeric or character chromosome numbers in the first column, numeric local probe positions in the second, and numeric copy number data for one or more samples in subsequent columns. The header of the copy number columns should be the sample IDs.}
  \item{segments}{a data frame or a list of data frames containing the segmentation results found by either \code{\link{pcf}} or \code{\link{multipcf}}.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{sample}{a numeric vector indicating which sample(s) is (are) to be plotted. The number(s) should correspond to the sample's place (in order of appearance) in \code{data}, or in \code{segments} in case \code{data} is unspecified.}
   \item{assembly}{a string specifying which genome assembly version should be applied to define the chromosome ideogram. Allowed options are "hg19", "hg18", "hg17" and "hg16" (corresponding to the four latest human genome annotations in the UCSC genome browser).}
   \item{winsoutliers}{an optional data frame of the same size as \code{data} identifying observations classified as outliers by \code{\link{winsorize}}. If specified, outliers will be marked by a different color and symbol than the other observations (see \code{wins.col} and \code{wins.pch}).}
  \item{xaxis}{either "pos" or "index". The former implies that the xaxis will represent the genomic positions, whereas the latter implies that the xaxis will 
  represent the probe index. Default is "pos".}
  \item{layout}{an integer vector of length two giving the number of rows and columns in the plot. Default is \code{c(1,1)}.}
  \item{\dots}{other graphical parameters. These include the common plot arguments \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}, \code{col} (default is "grey"), \code{pch} (default is 46, equivalent to "."), \code{cex}, \code{cex.lab}, \code{cex.main}, \code{cex.axis}, \code{las}, \code{tcl}, \code{mar} and \code{mgp} (see \code{\link{par}}
  on these). In addition, a range of graphical arguments specific for copy number plots may be specified, see \code{\link{plotSample}} on these.}
}

\details{
Several plots may be produced on the same page with the \code{layout} option. If the number of plots exceeds the desired page layout, the user is prompted before advancing to the next page of output.  
}
\note{
This function applies \code{par(fig)}, and is therefore not compatible with other setups for arranging multiple plots in one device such as \code{par(mfrow,mfcol)}. 
}

\seealso{
\code{\link{plotSample}}, \code{\link{plotChrom}}
}


\author{Gro Nilsen}
\examples{
#Lymphoma data
data(lymphoma)
#Take out a smaller subset of 6 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:6)

#Winsorize data:
wins.data <- winsorize(data=sub.lymphoma,return.outliers=TRUE)

#Use pcf to find segments:        
uni.segments <- pcf(data=wins.data,gamma=12)

#Use multipcf to find segments as well:
multi.segments <- multipcf(data=wins.data,gamma=12)

#Plot data and pcf-segments over entire genome for all six samples (one page
#for each sample):
plotGenome(data=sub.lymphoma,segments=uni.segments)

#Let each sample define its own range, and adjust range to fit all observations:
plotGenome(data=sub.lymphoma,segments=uni.segments,equalRange=FALSE,q=0)

#Add results from multipcf on top for four of the samples and let all plots 
#show on one page:
plotGenome(data=sub.lymphoma,segments=list(uni.segments,multi.segments),
    layout=c(2,2),sample=c(1:4))
    
#Change segment-colors, line widths, and legend:
plotGenome(data=sub.lymphoma,segments=list(uni.segments,multi.segments),layout=c(2,2),
    seg.col=c("red","blue"),seg.lwd=c(3,2),legend=c("uni","multi")
    ,sample=c(1:4))
    
#Aberration calling may be done by defining thresholds that determines the cuf-off 
#for what should be considered biologically significant aberrations. In this 
#example segments which are above 0.2 or below -0.2 are considered aberrated
#regions:     
plotGenome(segments=uni.segments,sample=5,connect=FALSE)
abline(h=0.2,col="blue",lty=5)
abline(h=-0.2,col="blue",lty=5)

}



