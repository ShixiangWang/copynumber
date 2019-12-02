
\name{subsetData}
\alias{subsetData}

\title{
Retrieve a data subset 
}
\description{
This function returns a subset of copy number data according to the input and the specified chromosomes and/or samples.
}
\usage{
 subsetData(data, chrom = NULL, sample = NULL, sep="\t", ...)
}
\arguments{
  \item{data}{either a data frame or the name of a tab-separated file from which copy number data can be read. The rows of the data frame or file should represent the probes. Column 1 must hold numeric or character chromosome numbers, column 2 the numeric local probe positions, and subsequent columns the numeric copy number measurements for one or more samples.}
  \item{chrom}{a numeric or character vector with chromosome(s) for which data should be selected. If unspecified, all chromosomes in \code{data} will be selected.}
  \item{sample}{a numeric vector indicating for which sample(s) data should be selected. The number(s) should correspond to the sample's place (in order of appearance) in \code{data}.}
  \item{sep}{the separator of the input files if \code{data}. Default is tab.}
  \item{\dots}{optional parameters to be passed to \code{read.table} in the case where \code{data} is to be read from file.}
}

\value{ 
A data frame containing the desired subset of data.

}
\author{Gro Nilsen}

\examples{
#Load lymphoma data
data(lymphoma)

#Select data only for samples 1 and 6 and chromosomes 1:9:
sub.data <- subsetData(data=lymphoma,chrom=c(1:9),sample=c(1,6))

}
