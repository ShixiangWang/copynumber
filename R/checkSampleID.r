
####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################


#Function to check that samples are found in both data and segments 
#If sample=NULL; all unique sampleIDs found in both data and segments are selected

#Input:
### data: the copy number data frame
### segments: a list containing one or more segmentations
### sample: a vector of sample numbers, or NULL


##Output:
### sampleID: a vector of unique selected sampleIDs to be plotted

##Function called by:
### checkAndRetrievePlotInput

##Function calls: none

checkSampleID <- function(data,segments,sample){
 
  sampleID <- NULL
  
  if(!is.null(data)){
    sample.names <- colnames(data)[-c(1:2)]
    #Pick out ID for samples to be plotted
    if(!is.null(sample)){
      sampleID <- sample.names[sample]
      if(length(sampleID)==0){
        stop("Input in 'sample' is larger than the number of samples found in 'data'",.call=FALSE)
      }
    }else{
      sampleID <- sample.names
    }
  }
  
  if(!is.null(segments)){
    all.segid <- sapply(segments,"[",i=1)   #returns a list with all sampleID for each segmentation
    all.segid <- lapply(all.segid,as.character)   #convert from factor to charachter
    all.segid <- lapply(all.segid,unique)         #get unique sampleids in each segmentation

    #Find common sampleIDs in data and all segmentations, given input in sample
    i <- 1
    if(is.null(sampleID)){    #data is NULL
      #Get sampleIDs to be plotted for the first segmentation result
      sampleID <- all.segid[[1]] 
      if(!is.null(sample)){
        sampleID <- sampleID[sample]
        if(length(sampleID)==0){
          stop("Input in 'sample' is larger than the number of samples found in 'segments'",.call=FALSE)
        }
      }
      i <- i+1
    } 
    #Get the sampleIDs that are also found in other segmentations (if more than one) 
    while(i<=length(all.segid)){
      sampleID <- intersect(sampleID,all.segid[[i]])
      i <- i+1
    }
    if(length(sampleID)==0){
      if(!is.null(data)){
        stop("no sampleIDs are common in 'data' and 'segments'",.call=FALSE)
      }else{
        stop("no sampleIDs are common in all components of 'segments'",.call=FALSE)
      }
    }
  }  
  
  
  #Check input sampleID and print errors or warnings if necessary:  
  sampleID <- sampleID[!is.na(sampleID)]  #could be NA if 'sample' is outside the number of samples represented in data/segments
      
  return(sampleID)

}#end checkSampleID

