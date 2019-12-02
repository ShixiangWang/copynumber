\name{plotAllele}                                      
\alias{plotAllele}

\title{Plot SNP data and/or aspcf segmentation results}

\description{
 Plot bivariate SNP data and/or aspcf segmentation results for each sample separately with chromosomes in different panels
}
\usage{
plotAllele(logR = NULL, BAF = NULL, segments = NULL, pos.unit = "bp", 
            sample = NULL, chrom = NULL, assembly="hg19", baf.thres = 
            c(0.1,0.9), winsoutliers = NULL, xaxis = "pos", layout = c(1,1), 
            plot.ideo = TRUE, ...)

}
\arguments{
  \item{logR}{a data frame with numeric or character chromosome numbers in the first column, numeric local probe positions in the second, and numeric copy number data for one or more samples in subsequent columns. The header of the copy number column(s) should give the sample IDss.}
  \item{BAF}{a data frame on the same format and size as \code{logR}, with chromosomes and local probe positions in the two first columns, and numeric BAF-measurements for one or more samples in subsequent columns.}
  \item{segments}{a data frame or a list of data frames containing the segmentation results found by \code{\link{aspcf}}.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{sample}{a numeric vector indicating which sample(s) is (are) to be plotted. The number(s) should correspond to the sample's place (in order of appearance) in \code{logR}, or in \code{segments} if \code{logR} is unspecified.}
  \item{chrom}{a numeric or character vector with chromosome number(s) to indicate which chromosome(s) is (are) to be plotted.}
  \item{assembly}{a string specifying which genome assembly version should be applied to define the chromosome ideogram. Allowed options are "hg19", "hg18", "hg17" and "hg16" (corresponding to the four latest human genome annotations in the UCSC genome browser).}
  \item{baf.thres}{a numeric vector of length 2 giving thresholds below/above which BAF-values will not be plotted (use this to remove germline homozygous BAF probes from the plot).}
  \item{winsoutliers}{an optional data frame of the same size as \code{logR} identifying observations classified as outliers by \code{\link{winsorize}}. If specified, outliers will be marked by a different color and symbol than the other observations (see \code{wins.col} and \code{wins.pch}).}
  \item{xaxis}{either "pos" or "index". The former implies that the xaxis will represent the genomic positions, whereas the latter implies that the xaxis will represent the probe indices. Default is "pos".}
  \item{layout}{an integer vector of length two giving the number of rows and columns in the plot. Default is \code{c(1,1)}.}
  \item{plot.ideo}{a logical value indicating whether the chromosome ideogram should be plotted. Only applicable when \code{xaxis="pos"}.}
  \item{\dots}{other graphical parameters. These include the common plot arguments \code{xlab}, \code{ylab}, \code{main}, \code{xlim}, \code{ylim}, \code{col} (default is "grey"), \code{pch} (default is 46, equivalent to "."), \code{cex}, \code{cex.lab}, \code{cex.main}, \code{cex.axis}, \code{las}, \code{tcl}, \code{mar} and \code{mgp} (see \code{\link{par}}
  on these). In addition, a range of graphical arguments specific for copy number plots may be specified, see \code{\link{plotSample}} on these. }
}

\details{
Several chromosome may be displayed on the same page with the \code{layout} option. If the number of chromosomes exceeds the desired page layout, the user is prompted before advancing to the next page of output.  
}


\note{
  This function applies \code{par(fig)}, and is therefore not compatible with other setups for arranging multiple plots in one device such as \code{par(mfrow,mfcol)}. 
}

\author{Gro Nilsen}
\examples{
#Load logR and BAF data:
data(logR)
data(BAF)

#Run aspcf::
aspcf.segments <- aspcf(logR,BAF)

#Plot
plotAllele(logR,BAF,aspcf.segments,layout=c(2,2))
}



