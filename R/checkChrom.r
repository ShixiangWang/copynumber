####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to check that specified chrom are numeric, found in data and/or segments and print appropriate error/warning messages depending on the plot type
#If chrom=NULL; all unique chromosomes found in data and/or segments are selected

#Input:
### chrom: a vector of unique chromosome numbers, or NULL
### data: the copy number data frame
### segments: a list containing one or more segmentations
### type: plot type (sample,chromosome,aspcf)


##Output:
### use.chrom : a vector of unique selected chromosomes to be plotted

##Function called by:
### plotFreq
### plotHeatmap
### plotWeightedFreq
### checkAndRetrievePlotInput

##Function calls:
### numericChrom

checkChrom <- function(data,segments,chrom){

  chrom.list <- NULL
  
  if(!is.null(data)){
    chrom.list <- unique(data[,1])
  }
  if(!is.null(segments)){
    all.chrom <- sapply(segments,"[",i=2)   #returns a list with chromsome numbers for each segmentation
    seg.chrom <- lapply(all.chrom,unique)
    
    
    #Find common chromosomes in data and all segmentations:
    i <- 1
    if(is.null(chrom.list)){
      chrom.list <- seg.chrom[[1]]
      i <- i+1
    }
    while(i<=length(seg.chrom)){
      chrom.list <- intersect(chrom.list,seg.chrom[[i]])
      i <- i+1
    }
  }#endif

  #Check input chromosomes, print error or warning if needed:
  if(is.null(chrom)){
      use.chrom <- chrom.list
  }else{
    #Make sure chrom is numeric
    chrom <- numericChrom(chrom)
    
    #Check that specified chrom are all in allowed list (found in data and segments):
    use <- chrom%in%chrom.list
    use.chrom <- chrom[use]
    if(length(use.chrom)==0){
      #print error message
      if(!is.null(data) && !is.null(segments)){
        msg <- "Specified chromosome(s) are not found in data and/or segments"
      }else if(is.null(segments)){
        msg <- "Specified chromosome(s) are not found in data"
      }else{
        msg <- "Specified chromosome(s) are not found in segments"
      }
      stop(msg,call.=FALSE)
    }else if(length(use.chrom)!=length(chrom)){
        not.use <- paste(chrom[!use],sep="",collapse=",")
        if(!is.null(data) && !is.null(segments)){
          msg <- "The following chromosome(s) are not found in data and/or segments:"
        }else if(is.null(segments)){
          msg <- "The following chromosome(s) are not found in data:"
        }else{
          msg <- "The following chromosome(s) are not found in segments:"
        }
        warning(paste(msg,not.use,sep=" "),call.=FALSE,immediate.=TRUE)
    }
  }#endif
  
  #make sure that chromosome numbers are in increasing order	
  use.chrom <- sort(use.chrom)       
  
  return(use.chrom)


}#end checkChrom



