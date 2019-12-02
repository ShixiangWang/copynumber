####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

### Function that pulls out the relevant contents of the input in res
### Makes it possible for res input to come directly from the segmentation algorithms or winsorize algorithm (regardless of whether res contains only a data frame or a list with segmentation/winsorize results).

##Required by:
### callAberrations
### imputeMissing
### interpolate.pcf
### checkAndRetrievePlotInput
### checkWinsoutliers
### checkSegments
### selectSegments
### subsetSegments
### plotCircle
### plotFreq
### plotHeatmap
### subsetData
### pcf
### pcfPlain
### multipcf
### aspcf


#what could be "segments","estimates","wins.data" or "wins.outliers"

pullOutContent <- function(res,what="segments"){

  #check if input is data frame or list
  if(!is.data.frame(res)){
    #res could either be a list containing the two segmentation elements segments and estimates, a list containing the two winsorize elements wins.data and wins.outliers, a list containing several segments as data frames, or a list containing several lists with segmentation results
    if("segments" %in% names(res)){
      #can assume that the list contains the output from one of the segmentation algorithms and has names segments and estimates
      #pick out the desired data frame depending on input in what:
      if(what=="estimates"){
        if("logR_estimates" %in% names(res)){
          #Segmentation results come from aspcf which returns a different name for estimates
          res <- res$logR_estimates  
        }else{
          res <- res$estimates
        }
      }else if(what=="segments"){
        res <- res$segments
      }
    }else if("wins.data" %in% names(res)){
      #can assume that the list contains output from winsorize and has names wins.data and wins.outliers
      #pick out the desired data frame depending on input in what:
      if(what %in% c("wins.data","estimates")){
        res <- res$wins.data
      }else if(what=="wins.outliers"){
        res <- res$wins.outliers
      }
    }else{
      #assume that the list contains different segmentation results
      #if each element in the list is itself a list containing segments and estimates need to pull out the segments (functions that take estimates as input does not have the option of specifying list of estimates)
      nSeg <- length(res)
      for(l in 1:nSeg){
        if(!is.data.frame(res[[l]])){
          #pull out the segments data frame
          res[[l]] <- res[[l]]$segments
        }    
      }
    } 
  }
  #If a data frame, res should be already be a data frame with the correct information, and does not need modification
   
  return(res)
  
}
