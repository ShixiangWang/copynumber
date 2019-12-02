\name{plotSample}                                      
\alias{plotSample}

\title{Plot copy number data and/or segmentation results by sample}
\description{
 Plot copy number data and/or segmentation results for each sample separately with chromosomes in different panels. 
}
\usage{
plotSample(data = NULL, segments = NULL, pos.unit = "bp", sample = NULL,
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
  on these). In addition, a range of graphical arguments specific for \code{plotSample} (as well as the similar functions \code{plotChrom}, \code{plotGenome} and \code{plotAllele}) may be specified: 
  	\describe{
  	  \item{\code{dir.print}:}{an optional directory where the plot(s) is (are) to be saved as pdf file(s). Defaults to NULL which implies that the plot(s) is (are) printed to screen instead.}
      \item{\code{file.name}}{an optional character vector containing file name(s) for the pdf file(s) to be saved.}
      \item{\code{onefile}:}{logical value indicating whether all plots should be plotted in one device / saved in one file. Default is TRUE. If FALSE, a new window is opened or a new file is saved for each sample (each chromosome for \code{plotChrom}).}
  	  \item{\code{plot.size}:}{a numeric vector of length 2 giving the width and height of the plotting window. Default is \code{c(11.6,8.2)}.}
      \item{\code{title}:}{an overall title for all plots on one page.}
      \item{\code{plot.unit}:}{the desired unit to be applied for probe position tick marks along the x-axis. Only "mbp" (default) and "kbp" is allowed.}
		  \item{\code{equalRange}:}{logical value indicating whether the range of the y-axis should be the same across all plots. Defaults to TRUE.}
		  \item{\code{q}:}{a numerical value in the range 0 to 1 indicating that \code{ylim} will be set to only include observations between the (1-q/2)- and the (q/2)-quantile. Observations that fall outside these quantiles are truncated to the limits of the plot, and are by default marked by a special symbol (see \code{q.pch}). Default is \code{q=0.01} when \code{data} is specified, and \code{q=0} otherwise.}
		  \item{\code{q.col}, \code{wins.col}:}{colors used to plot truncated observations and outliers. Default is "grey" and "magenta", respectively.}
      \item{\code{q.pch}, \code{wins.pch}:}{symbols used to plot truncated observations and outliers. Default is 42 (equivalent to "*") for both. Note that input must be of the same class as \code{pch} (numeric or character).}
		  \item{\code{q.cex}, \code{wins.cex}:}{magnification used for truncated observations and outliers relative to \code{cex}. Default is 0.4 for both.}
		  \item{\code{h}:}{a numerical value indicating that a horizontal reference line should be plotted at \code{y=h}. Default is \code{h=0}. \code{h=NULL} suppresses
      the plotting of a reference line.}
       \item{\code{at.x}:}{the points at which tick-marks on x-axis are to be drawn.}
		  \item{\code{at.y}:}{the points at which tick-marks on y-axis are to be drawn.}
      \item{\code{main.line}:}{the margin line for the main title.}
		  %\item{\code{h.col}:}{colour of the reference line. Default is "darkgrey"}
		  %\item{\code{h.lty}:}{line type used for the reference line. Default is 5 (longdash)}
		  %\item{\code{h.lwd}:}{line width for the reference line}
		  \item{\code{legend}:}{either a logical value indicating whether legends should be added to the plot if there is more than one segmentation result present in \code{segments}, or a character vector giving the legend texts to be used for the segmentation results. Default is TRUE, in which case the legend will be plotted in the topright corner of each plot.}
		  \item{\code{seg.col}:}{color(s) used to plot the segmentation result(s). The default colors are found using the function \code{rainbow(n)}, where \code{n} is the number of segmentation results found in \code{segments} (see \code{\link{rainbow}} for details).}
		  \item{\code{seg.lty}:}{the line type(s) used to plot the segmentation result(s). Default is 1.}
		  \item{\code{seg.lwd}:}{the line width(s) used to plot the segmentation result(s).}
		  \item{\code{connect}:}{logical value indicating whether segments should be connected by vertical lines, default is TRUE.}
		  \item{\code{ideo.frac}:}{a numerical value in the range 0 to 1 indicating the fraction of the plot to be occupied by the chromosome ideogram.}
      \item{\code{cyto.text}:}{a logical value indicating whether cytoband-names should be plotted along with the ideogram. Not recommended when many plots are
      plotted in the same grid, default is FALSE.}
      \item{\code{cex.cytotext}:}{the magnification used for the plotting of the cytoband-names.}
      \item{\code{cex.chrom}:}{the text size used to plot chromosome numbers in \code{plotGenome}.}
	 }
  }
}

\details{
Several plots may be produced on the same page with the \code{layout} option. If the number of plots exceeds the desired page layout, the user is prompted before advancing to the next page of output.  
}
\note{
These functions apply \code{par(fig)}, and are therefore not compatible with other setups for arranging multiple plots in one device such as \code{par(mfrow,mfcol)}. 
}


\seealso{
\code{\link{plotChrom}}, \code{\link{plotGenome}}
}

\author{Gro Nilsen}
\examples{
#Lymphoma data
data(lymphoma)
#Take out a smaller subset of 6 samples (using subsetData):
sub.lymphoma <- subsetData(lymphoma,sample=1:6)

#Winsorize data:
wins.data <- winsorize(data=sub.lymphoma)

#Use pcf to find segments:        
uni.segments <- pcf(data=wins.data,gamma=12)

#Use multipcf to find segments as well:
multi.segments <- multipcf(data=wins.data,gamma=12)

#Plot data and pcf-segments for one sample separately for each chromosome:
plotSample(data=sub.lymphoma,segments=uni.segments,sample=1,layout=c(5,5))
#Add cytoband text to ideogram (one page per chromosome to ensure sufficient 
#space)
plotSample(data=sub.lymphoma,segments=uni.segments,sample=1,layout=c(1,1),
    cyto.text=TRUE)
#Add multipcf-segmentation results, drop legend
plotSample(data=sub.lymphoma,segments=list(uni.segments,multi.segments),sample=1,
    layout=c(5,5),seg.col=c("red","blue"),seg.lwd=c(3,2),legend=FALSE)
#Plot by chromosome for two samples, but only chromosome 1-9. One window per 
#sample:
plotSample(data=sub.lymphoma,segments=list(uni.segments,multi.segments),sample=
    c(2,3),chrom=c(1:9),layout=c(3,3),seg.col=c("red","blue"),
    seg.lwd=c(3,2),onefile=FALSE)

#Zoom in on a particular region by setting xlim:    
plotSample(data=sub.lymphoma,segments=uni.segments,sample=1,chrom=1,plot.ideo=
    FALSE,xlim=c(140,170))
    
}



