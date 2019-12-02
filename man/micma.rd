\name{micma}
\alias{micma}

\docType{data}
\title{
Subset of 244K aCGH data
}
\description{
A subset of the 244K MicMa data set containing copy number measurements for six samples on chromosome 17. 
}
\usage{
data(micma)
}
\format{
 Data frame containing 7658 probes with log2-ratio copy numbers for 6 samples on chromosome 17. The first column contains the chromosome numbers, the second gives the local probe positions (in base pairs), while the subsequent columns contain the copy number measurements for the individual samples.
}
\source{
Mathiesen et al., "High resolution analysis of copy number changes in disseminated tumor cells of patients with breast cancer", Int J Cancer 131(4):E405:E415, 2011
}

\examples{
#Get data
data(micma)
}
