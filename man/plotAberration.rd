\name{plotAberration}
\alias{plotAberration}
\title{
Plot areas with copy number aberrations
}
\description{
Create plots reflecting the location of aberrated segments. Results may be visualized over the entire genome or by chromosomes. 
}
\usage{
plotAberration(segments, thres.gain, thres.loss = -thres.gain, pos.unit = "bp", 
    chrom = NULL, layout = c(1, 1),...)
}
\arguments{
  \item{segments}{a data frame containing the segmentation results found by either \code{\link{pcf}} or \code{\link{multipcf}}.}
  \item{thres.gain}{a numeric vector giving the threshold value(s) to be applied for calling gains.}
  \item{thres.loss}{a numeric vector of same length as \code{thres.gain} giving the threshold value(s) to be applied for calling losses. Default is to use the negative value of \code{thres.gain}.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{chrom}{a numeric or character vector with chromosome number(s) to indicate which chromosome(s) is (are) to be plotted. If unspecified the whole genome is plotted.}
  \item{layout}{the vector of length two giving the number of rows and columns in the plot window. Default is \code{c(1,1)}.}
  \item{\dots}{other optional graphical parameters. These include the plot arguments \code{xlab}, \code{ylab}, \code{main}, 
  \code{cex.main}, \code{mgp}, \code{cex.lab}, \code{cex.axis}, \code{mar} and \code{title} (see \code{\link{par}} on these), as well as \code{plot.size}, \code{plot.unit}, \code{plot.ideo}, \code{ideo.frac}, \code{cyto.text}, \code{assembly} and \code{cex.cytotext} (see \code{\link{plotSample}} on these). In addition, a range of graphical arguments 
  specific for this plot function may be specified: 
  	\describe{
  		\item{\code{colors}}{a character vector of length two where the first and second element specifies the color used to represent loss and gain, respectively. Default is c("dodgerblue","red").}
  		\item{\code{sample.labels}}{a logical value indicating whether sample labels are to be plotted along the y-axis. Default is TRUE.}
  		\item{\code{sep.samples}}{a number in the range 0 to 0.4 used to create some space between samples. 0 implies that there is no space. Default is 2/nsample, where nsample is the number of samples found in \code{segments}.}
  		\item{\code{sample.line}}{a numeric scalar giving the margin line where the sample labels should be written, starting at 0 counting outwards. Default is 0.2.}
  		\item{\code{sample.cex}}{the size of the sample labels.}
  	}
  }
}
\details{
For each sample, the aberrated regions are shown in the color specified in \code{colors[1]} (default dodgerblue) if the segment value is below \code{thres.loss} and the color specified in \code{colors[2]} (default red) if the segment value is above \code{thres.gain}. Non-aberrated regions are shown in white. Each row in the plot represents a sample, while probe positions are reflected along the x-axis.
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

#Plot aberrations for the entire genome
plotAberration(segments=seg,thres.gain=0.15)

#Plot aberrations for the first 4 chromosomes:
plotAberration(segments=seg,thres.gain=0.1,chrom=c(1:4),layout=c(2,2))

}
