
\name{subsetSegments}
\alias{subsetSegments}

\title{
Retrieve a subset of segments
}
\description{
This function returns a subset of segments according to the input and the specified chromosomes and/or samples.
}
\usage{
 subsetSegments(segments, chrom = NULL, sample = NULL, sep="\t", ...)
}
\arguments{
  \item{segments}{either a data frame or the name of a tab-separated file from which segmentation results can be read. Segmentation results may come from \code{\link{pcf}}, \code{\link{multipcf}} or \code{\link{aspcf}}.}
  \item{chrom}{a numeric or character vector with chromosome(s) for which segments should be selected. If unspecified, all chromosomes in \code{segments} will be selected.}
  \item{sample}{a numeric vector indicating for which sample(s) segments should be selected. The number(s) should correspond to the sample's place (in order of appearance) in \code{segments}.}
  \item{sep}{the separator of the input files if \code{segments} is a file. Default is tab.}
  \item{\dots}{optional parameters to be passed to \code{read.table} in the case where \code{segments} is to be read from file.}
}

\value{ 
A data frame containing the desired subset of segments. 
}
\author{Gro Nilsen}

\examples{
#Load lymphoma data
data(lymphoma)

#Select segments only for samples 1 and 6 and chromosomes 1:9:
segments <- pcf(lymphoma,gamma=12)
sub.segments <- subsetSegments(segments=segments,chrom=c(1:9),sample=c(1,6))

}
