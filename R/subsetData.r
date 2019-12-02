####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Get a subset of data for given chromosomes and/or sampleIDs

##Input:
### data : either a numeric matrix/data frame or the name of a file which data can be read from. The matrix or file should have chromosome numbers in the first column, local probe positions in the second, and copy number data for one or more samples in subsequent columns, and the header of copy number column(s) should be a sample identifier
### chrom : a numeric vector with chromosome number(s) for which data or segments should be selected. If unspecified, all chromosomes in data or segments will be selected
### sample : a numeric vector indicating which sample(s) is to selected. The number(s) should correspond to the sample's place (in order of appearance) in the data 
### sep : the separator of the input files if \code{data} is a file. Default is "\t"
### ... : optional parameters to be passed to \code{read.table} in the case where data are to be read from files

## Output:
### sel.data : subset of data containing required chromosomes and/or samples

##Required by:
### plotGamma

##Requires:
### numericChrom
### pullOutContent


subsetData <- function(data,chrom=NULL,sample=NULL,sep="\t",...){

  #Check if data is a file:
  isfile <- class(data)=="character"
  
  #get header and chrom from data
  if(isfile){
    head <- scan(data,nlines=1,what="character",quiet=TRUE,sep=sep) #Read header
    #Read just the two first columns to get chrom and pos
    data.chrom <- read.table(file=data,sep=sep,header=TRUE,colClasses=c(NA,rep("NULL",length(head)-1)),as.is=TRUE)[,1]   #could be character or numeric
  }else{
    #In case data comes from winsorize: check and possibly pull out data frame with wins.data
    data <- pullOutContent(data,what="wins.data")
    head <- colnames(data)
    data.chrom <- data[,1]
  }
  
  #check
  if(length(head)<3){
    stop("Data must have at least 3 columns",call.=FALSE)
  }
  data.sampleid <- head[-c(1:2)]
  nSample <- length(data.sampleid)
  #make sure chrom in data are numeric (covert x/X to 23 and y/Y to 24)
  data.chrom <- numericChrom(data.chrom)

  #Pick out relevant sampleIDs:
  keepsample <- 1:nSample
	if(!is.null(sample)){
    keepsample <- keepsample[sample]
    if(any(is.na(keepsample))){
      stop("Input in 'sample' is outside the number of samples represented in data",.call=FALSE)
    }
	}
 
	#Pick out relevant chromosomes:

	keepchrom <- rep(TRUE,length(data.chrom))
	if(!is.null(chrom)){
    #make sure selected chrom are numeric:
    chrom <- numericChrom(chrom)
    #Check that these chrom are found in data.chrom
    use <- chrom%in%unique(data.chrom)
    use.chrom <- chrom[use]
    if(length(use.chrom)==0){
      msg <- "Specified chromosome(s) is not found in data"
      stop(msg,call.=FALSE)
    }else if(length(use.chrom)!=length(chrom)){
      not.use <- paste(chrom[!use],sep="",collapse=",")
      msg <- paste("The following chromosome(s) are not found in data:",not.use,sep=" ")
      warning(msg,call.=FALSE,immediate.=TRUE)
    }
    keepchrom <- data.chrom%in%use.chrom
	}
	
	if(isfile){
    #Find the number of lines in file to read and to skip:
    d <- diff(keepchrom)
    cp <- which(d!=0)
    skip <- c()
    if(keepchrom[1]){
      #know that we should start reading at line 1
      skip <- c(skip,0)
    }
    skip <- c(skip,cp[d[cp]==1])

    nrows <- cp[d[cp]==-1]
    if(keepchrom[length(keepchrom)]){
      #know that we should end reading at line n
      nrows <- c(nrows,length(keepchrom))
    }
    nrows <- nrows - skip

    #Read relevant data
    sel.data <- data.frame(matrix(nrow=0,ncol=length(keepsample)+2))
    cc <- rep("NULL",nSample+2)
    cc[c(1:2,keepsample+2)] <- NA    #decide on which columns to read
    nreads <- length(nrows)
    for(i in 1:nreads){
      sel.data <- rbind(sel.data,read.table(data,header=FALSE,sep=sep,skip=skip[i]+1,nrows=nrows[i],colClasses=cc,as.is=TRUE))
    }
  }else{
    sel.data <- data[keepchrom,c(1:2,keepsample+2)]
  }
  colnames(sel.data) <- c(head[c(1:2)],data.sampleid[keepsample])
  
  return(sel.data)

}

