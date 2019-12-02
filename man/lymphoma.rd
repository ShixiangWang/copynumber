\name{lymphoma}
\alias{lymphoma}
\docType{data}
\title{
3K aCGH data
}
\description{
A subset of the aCGH data set taken from the reference below.
}
\usage{
data(lymphoma)
}
\format{
 Data frame containing 3091 probes with log2-ratio copy numbers for 21 samples. The first column contains the chromosome numbers, the second gives the local probe positions (in base pairs), while the subsequent columns contain the copy number measurements for the individual samples.
}
\source{
Eide et al., "Genomic alterations reveal potential for higher grade transformation in follicular lymphoma and confirm parallel evolution of tumor cell clones", Blood 116:1489-1497, 2010
}

\examples{
#Get data
data(lymphoma)
}
