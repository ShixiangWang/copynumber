####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

### Function that checks and modfies the input provided in segments in plot functions

##Input:
### segments: data frame or list with data frames containing segmentation results
### type: plot type

##Output:
### new.segments: Checked and modified list of segments

##Required by:
### checkAndRetrievePlotInput


##Requires:
### getUnisegFormat
### is.multiseg
### numericChrom
### pullOutContent


checkSegments <- function(segments,type){

    #Check and pull out relevant information in segments
    segments <- pullOutContent(res=segments,what="segments")
    
    #Get on list form if segments is matrix or data frame:
    if(class(segments)!="list"){
			segments <- list(seg=segments)
		}

		nSeg <- length(segments)
    seg.names <- names(segments)
		if(is.null(seg.names)){
			seg.names <- paste("Seg",1:nSeg,"")
		}
		#Segments could be on uniseg or multiseg form, convert all to uniseg form:
		new.segments <- vector("list",nSeg)  #empty list
		for(s in 1:nSeg){
			use.seg <- segments[[s]]
			multi <- is.multiseg(use.seg)
			if(multi){
				#stopifnot(ncol(use.seg)>=6)
				if(ncol(use.seg)<6){
				  stop("segments format is not correct",call.=FALSE)
				}
        #Convert to uniseg format:
				use.seg <- getUnisegFormat(use.seg)
			}else{
			  if(type!="aspcf"){
		      if(ncol(use.seg)<7){
				    stop("segments format is not correct",call.=FALSE)
            #stopifnot(ncol(use.seg)>=7)
          }
        }else{
          if(ncol(use.seg)<8){
            stop("segments format is not correct",call.=FALSE)
          }   
        } 
			}
			#Make sure chromosomes in column 2 are numeric:
			use.seg[,2] <- numericChrom(use.seg[,2])
			#Make sure arms in column 3 are character:
			use.seg[,3] <- as.character(use.seg[,3])
			
			new.segments[[s]] <- use.seg
		}

		names(new.segments) <- seg.names

		return(new.segments)
		
}#end checkSegments


