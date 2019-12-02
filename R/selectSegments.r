####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# FUNCTION THAT RETURNS A SELECTION OF MULTIPCF SEGMENTS BASED ON A DESIRED CHARACHTERISTIC

##Input:
### segments: data frame with segmentation results from multipcf
### what: selection criteria; one of "variance","length" and "aberration"
### thres: optional threshold for selcting, default is NULL
### nseg: desired number of segments, only used if thres=NULL
### p: minimum proportion when what="aberration"
### large: logical - select segments where variance, length or mean value is large (TRUE) or small (FALSE) for what equal to variance, length and aberration respectively.

##Required by:
###  none

##Requires:
### is.multiseg
### pullOutContent


selectSegments <- function(segments,what="variance",thres=NULL,nseg=10,large=TRUE,p=0.1){

  #Check input
  if(! what %in% c("variance","length","aberration")){
    stop("'what' must be one of variance, length and aberration")
  }
  #Make sure segments is a data frame
  segments <- pullOutContent(res=segments,what="segments")
  
  if(!is.multiseg(segments)){
    stop("'segments' must be on the format resulting from running multipcf!")
  }
  
  if(is.null(thres) && nseg > nrow(segments)){
    nseg <- nrow(segments)
    warning("'nseg' is larger than number of rows in 'segments'. Returning all segments.",call.=FALSE) 
    return(segments)
  }else{
    sel.res <- switch(what,
                      variance = subset.var(segments,nseg,thres,large),
                      length = subset.length(segments,nseg,thres,large),
                      aberration = subset.abe(segments,nseg,thres,p,large))
  
    #Sort sel.seg according to chromosome numbers:
    sel.res$sel.seg <- sel.res$sel.seg[order(sel.res$sel.seg[,1]),]
  
    return(sel.res)
  }

}

subset.var <- function(segments,nseg,thres,large){
  
  #calculate variance across samples for each segment:
  seg.var <- apply(segments[,-c(1:5)],1,var)
  
  if(!is.null(thres)){
    #Find the segments with variance above thres
    if(large){
      sel.seg <- segments[seg.var > thres,]
      if(nrow(sel.seg)==0){
        warning(paste("none of the segments have variance above ",thres,". Returning empty data frame.",sep=""),call.=FALSE)
      }
    }else{
      sel.seg <- segments[seg.var <thres,]
      if(nrow(sel.seg)==0){
        warning(paste("none of the segments have variance below ",thres,". Returning empty data frame.",sep=""),call.=FALSE)
      }
    }
  }else{
    #Find the nseg segments with the highest variance
    if(large){
      sel.seg <- segments[order(seg.var,decreasing=TRUE)[1:nseg],]
    }else{
      sel.seg <- segments[order(seg.var,decreasing=FALSE)[1:nseg],]
    }
  } 

  return(list(sel.seg=sel.seg,seg.var=seg.var))
}

subset.length <- function(segments,nseg,thres,large){
  #Find length of each segment:
  L <- segments[,4] - segments[,3] + 1
  if(!is.null(thres)){
    if(large){
      #Pick out long segments:
      sel.seg <- segments[L > thres,]
      if(nrow(sel.seg)==0){
        warning(paste("none of the segments are longer than ",thres,". Returning empty data frame.",sep=""),call.=FALSE)
      }
    }else{
      #Pick out short segments:
      sel.seg <- segments[L < thres,]
      if(nrow(sel.seg)==0){
        warning(paste("none of the segments are shorter than ",thres,". Returning empty data frame.",sep=""),call.=FALSE)
      }
    }
  }else{
    if(large){
      sel.seg <- segments[order(L,decreasing=TRUE)[1:nseg],]
    }else{
      sel.seg <- segments[order(L,decreasing=FALSE)[1:nseg],]
    }
  }

  return(list(sel.seg=sel.seg,seg.length=L))
}

subset.abe <- function(segments,nseg,thres,p,large){
  
  if(!is.null(thres)){
    if(large){
      prop.ab <- rowMeans(segments[,-c(1:5)] > thres)    
    }else{
      prop.ab <- rowMeans(segments[,-c(1:5)] < thres)                             
    }
    
    sel.seg <- segments[prop.ab >= p,]
    if(nrow(sel.seg)==0){
      if(large){
        warning(paste("none of the segments have mean value above ",thres,"for minimum ",p*100,"% of the samples. Returning empty data frame.",sep=""),call.=FALSE)
      }else{
        warning(paste("none of the segments have mean value below ",thres,"for minimum ",p*100,"% of the samples. Returning empty data frame.",sep=""),call.=FALSE)
      }
    }
    return(list(sel.seg=sel.seg,seg.ab.prop=prop.ab))
  }else{
    if(large){
      q <- apply(segments[,-c(1:5)],1,quantile,probs=1-p,type=1)
      q.ord <- order(q,decreasing=TRUE)
    }else{
      q <- apply(segments[,-c(1:5)],1,quantile,probs=p,type=1)
      q.ord <- order(q,decreasing=FALSE)   
    }
    sel.seg <- segments[q.ord[1:nseg],] 
  
    return(list(sel.seg=sel.seg,seg.quantile=q))
  }  
}