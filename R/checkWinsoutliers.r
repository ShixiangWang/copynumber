####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

## Function that check that input in winsoutliers has same dimension as data (in plots)

##Input:
### winsoutliers: data frame with outliers statueses
### data: data frame with logR data


##Required by:
### checkAndRetrievePlotInput


##Requires:
### pullOutContent

checkWinsoutliers <- function(winsoutliers,data){
	
	  winsoutliers <- pullOutContent(winsoutliers,what="wins.outliers")
		#Check that winsoutliers has same dimension as data:
		if(!all(dim(winsoutliers)==dim(data))){
			stop("winsoutliers must have the same number of rows and columns as data",call.=FALSE)
		}
		return(winsoutliers)
}#end checkWinsoutliers
