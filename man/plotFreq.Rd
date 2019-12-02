\name{plotFreq}
\alias{plotFreq}

\title{
Plot percentage of samples with an aberration at a genomic position
}
\description{
Plot the percentage of samples that have an amplification or deletion at a genomic position. Amplifications/deletions correspond to copy number values that are above/below a pre-defined threshold. Frequencies may be plotted over the entire genome or separately for each chromosome.
}
\usage{
plotFreq(segments, thres.gain, thres.loss = -thres.gain, pos.unit = "bp", 
    chrom = NULL, layout = c(1, 1),...)
}
\arguments{
  \item{segments}{a data frame containing the segmentation results found by either \code{\link{pcf}} or \code{\link{multipcf}}.}
  \item{thres.gain}{a numeric vector giving the threshold value(s) to be applied for calling gains.}
  \item{thres.loss}{a numeric vector of same length as \code{thres.gain} giving the threshold value(s) to be applied for calling losses. Default is to use the negative value of \code{thres.gain}.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{chrom}{a numeric or character vector with chromosome number(s) to indicate which chromosome(s) is (are) to be plotted. If unspecified the whole genome is plotted, otherwise each specified chromosome is plotted in a separate panel.}
  \item{layout}{an integer vector of length two giving the number of rows and columns in the plot. Default is \code{c(1,1)}.}
  \item{\dots}{other graphical parameters. These include the common plot arguments \code{xlab}, \code{ylab},  \code{title}, \code{main}, 
      \code{cex.main}, \code{mgp}, \code{cex.lab}, \code{cex.axis}, \code{ylim}, \code{xlim}, and \code{las}  (see \code{\link{par}} on these), as well as \code{plot.size}, \code{plot.unit}, \code{plot.ideo}, \code{ideo.frac}, \code{cyto.text}, \code{assembly} and \code{cex.cytotext} (see \code{\link{plotSample}} on these). In addition, some other graphical arguments specific for this plot function may be specified: 
      	\describe{
      	 	\item{\code{col.gain}}{the color to be used for the gain frequencies, default is "red".}
      		\item{\code{col.loss}}{the color to be used for the loss frequencies, default is "blue".}
      		\item{\code{continuous}}{a logical value indicating whether the probe frequencies should be presented as continuous, i.e. the plotted probe frequency will start
          and end halfway between adjacent probes (except across arms when chromosomes are plotted and except across chromosomes when genome is plotted). Default is TRUE.}
          \item{\code{percentLines}}{either a logical value indicating if horizontal percentages lines should be plotted, or a numeric vector with percentages at which 
          such lines should be plotted. Default is TRUE.}
      	}
  }
}
\details{
The percentage of samples with an aberration is calculated and plotted for all genomic positions. Regions with gain or loss will be those where copy number values are above or below the values given in \code{thres.gain} and \code{thres.loss}, respectively.
}


\note{
  This function applies \code{par(fig)}, and is therefore not compatible with other setups for arranging multiple plots in one device such as \code{par(mfrow,mfcol)}. 
}

\author{
Gro Nilsen
}

\examples{
#load lymphoma data
data(lymphoma)
#Run pcf
seg <- pcf(data=lymphoma,gamma=12)

#Plot over entire genome, gain and loss thresholds are 0.1 and -0.1:
plotFreq(segments=seg,thres.gain=0.1)

#Plot by chromosomes, two sets of thresholds:
plotFreq(segments=seg,thres.gain=c(0.1,0.2), thres.loss=c(-0.05,-0.1), chrom=c(1:23),
layout=c(5,5))
  

}

