\name{plotHeatmap}
\alias{plotHeatmap}
\title{
Plot copy number heatmap
}
\description{
Heatmap reflecting the magnitude of estimated copy numbers relative to some pre-defined limits. Estimates may be obtained using \code{pcf} or \code{multipcf}, and results may be visualized over the entire genome or by chromosomes. 
}
\usage{
plotHeatmap(segments, upper.lim, lower.lim = -upper.lim, pos.unit = "bp", 
    chrom = NULL, layout = c(1, 1),...)
}
\arguments{
  \item{segments}{a data frame containing the segmentation results found by either \code{\link{pcf}} or \code{\link{multipcf}}.}
  \item{upper.lim}{a positive numeric vector giving the upper limits(s) to be applied. The colors in the heatmap will reflect the magnitude of the estimated copy numbers relative to this limit, see details.}
  \item{lower.lim}{a negative numeric vector of same length as \code{upper.lim} giving the lower limits(s) to be applied. Default is to use the negative value of \code{upper.lim}.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{chrom}{a numeric or character vector with chromosome number(s) to indicate which chromosome(s) is (are) to be plotted. If unspecified the whole genome is plotted.}
  \item{layout}{the vector of length two giving the number of rows and columns in the plot window. Default is \code{c(1,1)}.}
  \item{\dots}{other optional graphical parameters. These include the plot arguments \code{xlab}, \code{ylab}, \code{main}, 
  \code{cex.main}, \code{mgp}, \code{cex.lab}, \code{cex.axis}, \code{mar} and \code{title} (see \code{\link{par}} on these), as well as \code{plot.size}, \code{plot.unit}, \code{plot.ideo}, \code{ideo.frac}, \code{cyto.text}, \code{assembly} and \code{cex.cytotext} (see \code{\link{plotSample}} on these). In addition, a range of graphical arguments 
  specific for this plot function may be specified: 
  	\describe{
  		\item{\code{colors}}{a character vector of length three giving the colors to interpolate in the heatmap, default is c("dodgerblue","black","red").}
  		\item{\code{n.col}}{an integer giving the number of color shades to be applied in the interpolation, default is 50.}
  		\item{\code{sample.labels}}{a logical value indicating whether sample labels are to be plotted along the y-axis. Default is TRUE.}
 		  \item{\code{sep.samples}}{a number in the range 0 to 0.4 used to create some space between samples. Default is 0, which implies that there is no space.}
  		\item{\code{sample.line}}{a numeric scalar giving the margin line where the sample labels should be written, starting at 0 counting outwards. Default is 0.2.}
  		\item{\code{sample.cex}}{the size of the sample labels.}
  	}
  }
}
\details{
For each sample, the segments are represented by a rectangle plotted in a color corresponding to the difference between the segment copy number value and the limits.
If the value is below \code{lower.lim}, the color of the rectangle will equal the input in \code{colors[1]} (default dodgerblue). If the value is above \code{lower.lim}, but below zero, the color of the rectangle will be a nuance between the input in \code{colors[1]} and \code{colors[2]} (default black). The closer the value is to zero, the closer the nuance will be to \code{colors[2]}. Similary, if the value is above \code{upper.lim}, the color of the rectangle will equal the input in \code{colors[3]} (default red), whereas if the value is below \code{upper.lim}, but above zero, the color will be a nuance between the input in \code{colors[2]} and \code{colors[3]}. Again, the closer the value is to zero, the closer the nuance will be to \code{colors[2]}.

Each row in the heatmap represents a sample, while probe positions are reflected along the x-axis.
}


\note{
  This function applies \code{par(fig)}, and is therefore not compatible with other setups for arranging multiple plots in one device such as \code{par(mfrow,mfcol)}. 
}

\author{
  Gro Nilsen
}
\examples{
#Load lymphoma data
data(lymphoma)

#Run pcf to obtain estimated copy number values
seg <- pcf(data=lymphoma,gamma=12)

#Heatmap for entire genome, two limit values:
plotHeatmap(segments=seg,upper.lim=c(0.1,0.5),layout=c(2,1))

#Heatmap for the first 4 chromosomes:
plotHeatmap(segments=seg,upper.lim=0.1,chrom=c(1:4),layout=c(2,2))

}
