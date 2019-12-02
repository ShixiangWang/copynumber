\name{SNPdata}
\alias{logR}
\alias{BAF}
\docType{data}
\title{
Artificial SNP array data
}
\description{
Artificial SNP array data containing a logR track and a BAF track
}
\usage{
data(logR)
data(BAF)
}
\format{
 Two corresponding data sets containing 10000 probes with logR and BAF measurements, respectively, for 2 samples. The two first columns in both data sets contain chromosome numbers and local probe positions (in base pairs), while the subsequent columns contain logR-values and BAF-values in the two data sets, respectively.
}


\examples{
#Get data
data(logR)
data(BAF)
}
