####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

###Function that checks plot input; and retrieves and returns necessary information

###Input:
### data: data.frame or matrix with observed copy numbers; chrom in col 1, pos in col 2 and copy numbers in remaining columns
### segments: segmentation results on either uniseg or multiseg format
### winsoutliers: outlier statuses from winsorization; same number of columns and rows as data
### type: type of plot; "sample","genome" or "chromosome"  
### xaxis: "index" eller "pos"
### pos.unit: "bp","kbp", "mbp"
### sample: the sample numbers to be plotted
### chrom: chromosomes to be plotted  (NB: if type=genome -> chrom=NULL)


### Output:
### data: checked data with numeric chromosomes in column 1
### segments: checked and modified segments according to checkSegments
### sampleID: a character vector with picked out sampleIDs to be plotted
### chrom: a vector with checked and picked out chromosomes to be plotted


##Required by:
### plotAllele
### plotChrom
### plotGenome
### plotSample

##Requires:
### checkChrom
### checkSegments
### checkWinsoutliers
### numericChrom
### pullOutContent

checkAndRetrievePlotInput <- function(data,segments,winsoutliers,type,xaxis,pos.unit,sample,chrom=NULL){

  if(type=="aspcf"){
    logR <- data$logR
    BAF <- data$BAF
    #Check that at least one of of logR, BAF and segments has been specified
    if(all(sapply(list(logR,BAF,segments),is.null))){
		  stop("Arguments 'logR' and 'BAF', or 'segments', must be specified!",call.=FALSE)  
    }
    if((is.null(logR)&&!is.null(BAF))||(is.null(BAF)&&!is.null(logR))){
      stop("'BAF' must be specified if 'logR' is specified, and vice versa",call.=FALSE)
    }
    if(!is.null(logR)){
      #logR and BAF could possibly have been winsorized: make sure it is a data frame:
      logR <- pullOutContent(logR,what="wins.data") 
      BAF <- pullOutContent(BAF,what="wins.data")
      #Check dimensions of logR and BAF:
      stopifnot(ncol(logR)>=3)
      stopifnot(ncol(logR)==ncol(BAF),nrow(logR)==nrow(BAF))
    }
    data <- logR    
  }else{
    #Check that at least of of data and segments has been specified
  	if(is.null(data)&&is.null(segments)){
  		stop("One of the arguments 'data' and 'segments' must be specified!",call.=FALSE)
  	}
	}
  	
	#Check and retrieve data input:
	if(!is.null(data)){
	  #data could possibly come from winsorize, make sure to pull out wins.data data frame
  	data <- pullOutContent(data,what="wins.data")
  	stopifnot(ncol(data)>=3)  #something is missing in input data
		#Make sure chromosomes are numeric:
    data[,1] <- numericChrom(data[,1])
    
    #Check winsoutliers
		if(!is.null(winsoutliers)){
      winsoutliers <- checkWinsoutliers(winsoutliers,data)
	  }

	}#endif

	#Check, modify and retrieve segments input:
	if(!is.null(segments)){
		segments <- checkSegments(segments,type)
	}#endif
	
	
	#Check sampleIDs to be plotted (only plot sampleIDs found in both data and segments)
	sampleID <- checkSampleID(data,segments,sample)

	#Check and if necessary modify chrom to be plotted:
	if(type!="genome"){
    chrom <- checkChrom(data,segments,chrom)
	}
	#Check xaxis input:
	if(!xaxis %in% c("pos","index")){
    stop("xaxis must be one of 'pos' and 'index'",call.=FALSE)
	}

	
	#Check that pos.unit has been specified if xaxis="pos"
	if(xaxis=="pos"){
    #Check pos.unit input:
    if(!pos.unit %in% c("bp","kbp","mbp")){
      stop("pos.unit must be one of bp, kbp and mbp",call.=FALSE)
    }
  }
  
  if(type=="aspcf"){
    #Return data as list
    data <- list(logR=data,BAF=BAF)
  }
  
  return(list(data=data,segments=segments,sampleID=sampleID,chrom=chrom,winsoutliers=winsoutliers))
} 

  
  
  
  
	
	