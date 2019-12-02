####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################


#Function that calls segments as gain, normal or loss

##Input:
### segments: segmentation results from pcf or multipcf
### thres.gain,thres.loss: threshold(s) to be applied for abberration calling

##Required by: none

##Requires: 
### is.multiseg
### pullOutContent


callAberrations <- function(segments,thres.gain,thres.loss=-thres.gain){
  
  #Make sure segments is a data frame
  segments <- pullOutContent(res=segments,what="segments")

  if(is.multiseg(segments)){
    call.seg <- matrix("normal",nrow=nrow(segments),ncol=ncol(segments)-5)
    colnames(call.seg) = colnames(segments)[-c(1:5)]
		gain <- segments[,-c(1:5),drop=FALSE] > thres.gain
		call.seg[gain] <- "gain"
    loss <- segments[,-c(1:5),drop=FALSE] < thres.loss
    call.seg[loss] <- "loss"
    
    call.seg <- as.data.frame(cbind(segments[,c(1:5)],call.seg))
    
  }else{
    call.seg <- rep("normal",nrow(segments))
		gain <- segments[,7] > thres.gain
		call.seg[gain] <- "gain"
    loss <- segments[,7] < thres.loss
    call.seg[loss] <- "loss"
    
    call.seg <- as.data.frame(cbind(segments[,c(1:6)],call.seg))
    colnames(call.seg)[7] = "call"
    
  }
  
  return(call.seg)

}

