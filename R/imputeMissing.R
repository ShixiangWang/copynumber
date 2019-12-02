####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# FUNCTION THAT IMPUTES MISSING VALUES IN COPY NUMBER DATA; EITHER A CONSTANT
# VALUE OR THE PCF-VALUE OF NEAREST NON-MISSING NEIGHBOUR

##Input:
### data : a numeric matrix with chromosome numbers in the first column, local probe positions in the second, and copy number data for one or more samples in subsequent columns
### method : the imputation method to be used. Must be one of "constant" and "pcf"
### c : a numerical value to be imputed if method is "constant". Default is 0
### pcf.est : a numeric matrix of the same size as data, with chromosome numbers and positions in the first two columns, and copy number estimates from \code{pcf} in the subsequent columns. Only applicable if \code{method="pcf"}. If unspecified, \code{pcf} is run internally to find estimates
### ... : other optional parameters to be passed to \code{pcf}

## Output:
### data matrix of same size as data with missing values imputed

##Required by:
### none

##Requires:
###  pcf
### pullOutContent



imputeMissing <- function(data,method,c=0,pcf.est=NULL,...){

  #Check input
  #First check if data comes from winsorize, and if so make sure it is a data frame
  data <- pullOutContent(data,what="wins.data")
  stopifnot(ncol(data)>=3)  #something is missing in input data
 
  #Check method input:
  if(!method %in% c("constant","pcf")){
    stop("method must be one of 'constant' and 'pcf'",call.=FALSE)
  }
  
	imp.data <- switch(method, 
				constant= imp.constant(data,c),
				pcf = imp.pcf(data,yhat=pcf.est,...))
				
	return(imp.data)
	
}#endimputeMissing


#Two possible imputation methods:
imp.constant <- function(data,c){
	cn.data <- data[,-c(1:2),drop=FALSE]
	na <- is.na(cn.data)
	imp.data <- cn.data
	imp.data[na] <- c
	
	imp.data <- cbind(data[,c(1:2),drop=FALSE],imp.data)
	return(imp.data)

}

imp.pcf <- function(data,yhat,...){
	
	na <- is.na(data)
	imp.data <- data
	
	if(is.null(yhat)){
		#Run pcf:
		s <- apply(na,2,any)  #find samples with missing values
		if(any(s[1:2]==TRUE)){
		  stop("missing values are not allowed in data columns 1 and 2", call.=FALSE)
		}
		pcf.res <- pcf(data=data[,c(1:2,which(s))],return.est=TRUE,...)  #Only run PCF on samples with missing values
		#Make sure yhat is the same size as cn.data to maintain location of na
		yhat <- data.frame(matrix(NA,nrow=nrow(data),ncol=ncol(data)))
		yhat[,c(1:2,which(s))] <- pcf.res$estimates
	}else{
    #Make sure pcf.est is a data frame
    yhat <- pullOutContent(res=yhat,what="estimates")
    if(ncol(yhat)!=ncol(data) || nrow(yhat)!=nrow(data)){
      stop("pcf.est must be the same size as data",call.=FALSE)
    }
	}
	imp.data[na] <- yhat[na]

	return(imp.data)
}



