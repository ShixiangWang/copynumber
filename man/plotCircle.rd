\name{plotCircle}
\alias{plotCircle}

\title{
Plot a circular genome with aberration frequencies and connections between genomic loci added. 
}
\description{
A circular genome is plotted and the percentage of samples that have a gain or a loss at a genomic position is added in the middle of the circle. Gains/losses correspond to copy number values that are above/below a pre-defined threshold. In addition arcs representing some connection between genomic loci may be added.
}
\usage{
plotCircle(segments, thres.gain, thres.loss = -thres.gain, pos.unit = "bp",
    freq.colors = c("red","limegreen"), alpha = 1/7, arcs = NULL, 
    arc.colors = c("goldenrod1","dodgerblue"), d = 0.3, assembly = "hg19")

}
\arguments{
  \item{segments}{a data frame containing the segmentation results found by either \code{\link{pcf}} or \code{\link{multipcf}}.}
  \item{thres.gain}{a scalar giving the threshold value to be applied for calling gains.}
  \item{thres.loss}{a scalar giving the threshold value to be applied for calling losses. Default is to use the negative value of \code{thres.gain}.}
  \item{pos.unit}{the unit used to represent the probe positions. Allowed options are "mbp" (mega base pairs), "kbp" (kilo base pairs) or "bp" (base pairs). By default assumed to be "bp".}
  \item{freq.colors}{a vector giving the colors to be used for the amplification and deletion frequencies, respectively. Default is c("red","limegreen").}
 	\item{alpha}{a scalar in the range 0 to 1 determining the amount of scaling of the aberration frequencies. For the default value of 1/7 the distance between the genome circle and the zero-line of the frequency-circle corresponds to an aberration percentage of 100 \%. See details. }
 	\item{arcs}{an optional matrix or data frame with 5 columns specifying connections between genomic loci. The first two columns must give the chromosome numbers and local positions for the start points of the arcs, while the two next columns give the chromosome numbers and local positions for the end point of arcs. The last column should contain a vector of numbers 1,2,... indicating that the arcs belong to different classes. Each class of arcs will then be plotted in a different color.}
 	\item{arc.colors}{a vector giving the colors to be used for plotting the different classes of arcs. Cannot be shorter than the number of unique classes in \code{arcs}. The first color will represent the first class in \code{arcs}, the second color the second class and so on.}
 \item{d}{a scalar > 0 representing the distance from the genome circle to the starting points of the arcs. Set d=0 to make arcs start at the genome circle.}
 \item{assembly}{a string specifying which genome assembly version should be applied to determine chromosome ideograms. Allowed options are "hg19", "hg18", "hg17" and "hg16" (corresponding to the four latest human genome annotations in the UCSC genome browser).}     		
}
\details{
 To zoom in on the observed aberration frequencies one may increase \code{alpha}. However, the user should be aware that this implies that the distance between the genome circle and the frequency zero-line does not reflect an aberration frequency of 100 \%. Since the distance between the two circles is always 1/7, the maximum plotted percentage will be 100/(alpha*7) and any percentages that are higher than this will be truncated to this value. 
}


\author{
Gro Nilsen
}

\examples{
#load lymphoma data
data(lymphoma)
#Run pcf
pcf.res <- pcf(data=lymphoma,gamma=12)

plotCircle(segments=pcf.res,thres.gain=0.1)

#Use alpha to view the frequencies in more detail:
plotCircle(segments=pcf.res,thres.gain=0.1,alpha=1/5)

#An example of how to specify arcs
#Using multipcf, we compute the correlation between all segments and then
#retrieve those that have absolute inter-chromosomal correlation > 0.7
multiseg <- multipcf(lymphoma)
nseg = nrow(multiseg)
cormat = cor(t(multiseg[,-c(1:5)]))
chr.from <- c()
pos.from <- c()
chr.to <- c()
pos.to <- c()
cl <- c()

thresh = 0.7
for (i in 1:(nseg-1)) {
  for (j in (i+1):nseg) {
    #Check if segment-correlation is larger than threshold and that the two 
    #segments are located on different chromosomes
    if (abs(cormat[i,j]) > thresh && multiseg$chrom[i] != multiseg$chrom[j]) {
      chr.from = c(chr.from,multiseg$chrom[i])
      chr.to = c(chr.to,multiseg$chrom[j])
      pos.from = c(pos.from,(multiseg$start.pos[i] + multiseg$end.pos[i])/2)
      pos.to = c(pos.to,(multiseg$start.pos[j] + multiseg$end.pos[j])/2)
      if(cormat[i,j] > thresh){
        cl <- c(cl,1)           #class 1 for those with positive correlation
      }else{
        cl <- c(cl,2)           #class 2 for those with negative correlation
      }    
    }
  }
}
  
arcs <- cbind(chr.from,pos.from,chr.to,pos.to,cl)  

#Plot arcs between segment with high correlations; positive correlation in
#orange, negative correlation in blue:
plotCircle(segments=pcf.res,thres.gain=0.15,arcs=arcs,d=0)


  

}

