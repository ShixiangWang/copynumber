####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that finds the y-range of segments to be plotted

#Input:
### segments: dataframe with a segmentation result
### equalRange: should range be found over all samples and chromosomes in segments?
### sampleID: a sampleID for which range is to be found
### chrom: a chromosome number for which range is to be found
### baf: are we looking for range of baf-segments? if so; these values are found in column 8, otherwise they are found in column 7

# Output:
### seg.lim: vector of length 2 giving the minimum and maximum segment limit

##Required by:
### plotAllele
### plotChrom
### plotGenome
### plotSample



##Requires:
### none

get.seglim <- function(segments,equalRange,sampleID=NULL,k=NULL,baf=FALSE){
	
	if(equalRange){
		#Use all segments to determine limits:
		use.segments <- segments
	}else{
		#Use only segments indicated by index or k to calculate limits:
		if(!is.null(k)){	
			use.segments <- segments[segments[,2]==k,]
		}else{if(!is.null(sampleID)){
			keep <- which(segments[,1]==sampleID)
			use.segments <- segments[keep,]
		}}
	}	
	seg.lim <- rep(NA,2)

  if(nrow(use.segments)>0){			
    if(!baf){
      seg.lim <- c(min(use.segments[,7]),max(use.segments[,7]))
    }else{
      seg.lim <- c(min(use.segments[,8]),max(use.segments[,8]))
    }
  }
	
	return(seg.lim)
	
}#end get.seglim