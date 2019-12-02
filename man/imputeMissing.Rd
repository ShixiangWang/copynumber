\name{imputeMissing}
\alias{imputeMissing}
\title{Impute missing copy number values}
\description{
  Missing copy number values are imputed by a constant value or pcf-estimates.
}
\usage{
imputeMissing(data, method, c = 0, pcf.est = NULL,...)
}
\arguments{
  \item{data}{a data frame with numeric or character chromosome numbers in the first column, numeric local probe positions in the second, and numeric copy number data for one or more samples in subsequent columns.}
  \item{method}{the imputation method to be used. Must be one of "constant" and "pcf".}
  \item{c}{a numerical value to be imputed if method is "constant". Default is 0.}
  \item{pcf.est}{a data frame of same size as \code{data}, with chromosome numbers and positions in the first two columns, and copy number estimates obtained from \code{pcf} in the subsequent columns. Only applicable if \code{method="pcf"}. If unspecified and \code{method="pcf"}, \code{pcf} is run internally to find estimates.}
  \item{\dots}{other relevant parameters to be passed on to \code{pcf}}
  }
\details{
The available imputation methods are:
\describe{
  \item{\code{constant}:}{all missing values in \code{data} are replaced by the specified value \code{c}.}
  \item{\code{pcf}:}{the estimates from pcf-segmentation (see \code{\link{pcf}}) are used to impute missing values. If \code{pcf} has already been run, these estimates may be specified in \code{pcf.est}. If \code{pcf.est} is unspecified, \code{pcf} is run on the input data. In \code{pcf} the analysis is done on the observed values, and estimates for missing observations are set to be the estimate of the nearest observed probe.}
  }
}
\value{
  A data frame of the same size and format as \code{data} with all missing values imputed.
}
\author{Gro Nilsen}

\seealso{\code{\link{pcf}}}

\examples{
#Load lymphoma data
data(lymphoma)
chrom <- lymphoma[,1]
pos <- lymphoma[,2]
#pick out data for the first six samples:
cn.data <- lymphoma[,3:8]

#Create missing values in cn.data at random positions:
n <- nrow(cn.data)*ncol(cn.data)
r <- matrix(rbinom(n=n,size=1,prob=0.95),nrow=nrow(cn.data),ncol=ncol(cn.data))
cn.data[r==0] <- NA    #matrix with approximately 5% missing values
mis.data <- data.frame(chrom,pos,cn.data)

#Impute missing values by constant, c=0:
imp.data <- imputeMissing(data=mis.data,method="constant")

#Impute missing values by obtained pcf-values:
pcf.est <- pcf(data=mis.data,return.est=TRUE)
imp.data <- imputeMissing(data=mis.data,method="pcf",pcf.est=pcf.est)

#Or run pcf within imputeMissing:
imp.data <- imputeMissing(data=mis.data,method="pcf")

}
