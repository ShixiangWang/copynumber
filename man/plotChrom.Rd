\name{plotChrom}                                      
\alias{plotChrom}

\title{Plot copy number data and/or segmentation results by chromosome}
\description{
 Plot copy number data and/or segmentation results for each chromosome separately with samples in different panels. 
}
\usage{
plotChrom(data = NULL, segments = NULL, pos.unit = "bp", sample = NULL,
            chrom = NULL, assembly = "hg19", winsoutliers = NULL, xaxis = 
            "pos", layout = c(1,1), plot.ideo = TRUE, ...)
}
\arguments{
  \item{data}{a data frame with numeric or character chromosome numbers in the first column, numeric local probe positions in the second, and numeric copy number data for one or more samples in subsequent columns. The header of the copy number columns should be the sample IDs.}
  \item{segments}{a data frame or a list of data frames containing the segmentation results found by either \code{\link{pcf}} or \code{\link{multipcf}}.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{sample}{a numeric vector indicating which sample(s) is (are) to be plotted. The number(s) should correspond to the sample's place (in order of appearance) in \code{data}, or in \code{segments} in case \code{data} is unspecified.}
  \item{chrom}{a numeric or character vector with chromosome number(s) to indicate which chromosome(s) is (are) to be plotted.}
  \item{assembly}{a string specifying which genome assembly version should be applied to define the chromosome ideogram. Allowed options are "hg19", "hg18", "hg17" and "hg16" (corresponding to the four latest human genome annotations in the UCSC genome browser).}
  \item{winsoutliers}{an optional data frame of the same size as \code{data} identifying observations classified as outliers by \code{\link{winsorize}}. If specified, outliers will be marked by a different color and symbol than the other observations (see \code{wins.col} and \code{wins.pch}).}
  \item{xaxis}{either "pos" or "index". The former implies that the xaxis will represent the genomic positions, whereas the latter implies that the xaxis will 
  represent the probe index. Default is "pos".}
  \item{layout}{an integer vector of length two giving the number of rows and columns in the plot. Default is \code{c(1,1)}.}
  \item{plot.ideo}{a logical value indicating whether the chromosome ideogram should be plotted. Only applicable when \code{xaxis="pos"}.}
  \item{\dots}{other graphical parameters. These include the common plot arguments \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}, \code{col} (default is "grey"), \code{pch} (default is 46, equivalent to "."), \code{cex}, \code{cex.lab}, \code{cex.main}, \code{cex.axis}, \code{las}, \code{tcl}, \code{mar} and \code{mgp} (see \code{\link{par}}
  on these). In addition, a range of graphical arguments specific for copy number plots may be specified, see \code{\link{plotSample}} on these. }
}

\details{
Several plots may be produced on the same page with the \code{layout} option. If the number of plots exceeds the desired page layout, the user is prompted before advancing to the next page of output.  
}
\note{
%In \code{plotChrom}, only chromosomes found in both \code{data} and all segmentation results in \code{segments} will be plotted.
This function applies \code{par(fig)}, and is therefore not compatible with other setups for arranging multiple plots in one device such as \code{par(mfrow,mfcol)}. 
}


\seealso{
\code{\link{plotSample}}, \code{\link{plotGenome}}
}

\author{Gro Nilsen}
\examples{
#Lymphoma data
data(lymphoma)
#Take out a smaller subset of 6 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:6)

#Winsorize data:
wins.res <- winsorize(data=sub.lymphoma,return.outliers=TRUE)

#Use pcf to find segments:        
uni.segments <- pcf(data=wins.res,gamma=12)

#Use multipcf to find segments as well:
multi.segments <- multipcf(data=wins.res,gamma=12)

#Plot data and segments for chromosome 1 separately for each sample:
plotChrom(data=sub.lymphoma,segments=list(uni.segments,multi.segments),chrom=1,
    layout=c(3,2))
#Let xaxis be probe index, and do not connect segments by vertical lines:
plotChrom(data=sub.lymphoma,segments=list(uni.segments,multi.segments),chrom=1,
    xaxis="index",layout=c(3,2),legend=FALSE,connect=FALSE)
#Data was winsorized earlier. Mark winsorized values by different color
#and symbol: 
plotChrom(data=wins.res,chrom=1,winsoutliers=wins.res,layout=c(3,2))
#Save plots in working directory:
\donttest{plotChrom(data=sub.lymphoma,segments=uni.segments,chrom=c(1,2),
    layout=c(3,2),dir.print=getwd(),file.name=c("chromosome1","chromosome2"),
    onefile=FALSE)}  
}



