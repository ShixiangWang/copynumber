####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

##Function to get multisegment results on a unisegmentformat:

##Input:
### segments: a data frame with segmentation results

##Output:
### uni.segments: the segments converted to unisegmentsformats; i.e sampleId in first column, chrom in second, arms in third, start and end pos in fourth and fifth, number of pos in sixth and mean value in seventh

##Required by:
### checkSegments

##Requires:
### is.multiseg


getUnisegFormat <- function(segments){

	#Check that the segments are really on a multiseg format first:
	stopifnot(is.multiseg(segments))
	
	nSample <- ncol(segments)-5
	nSeg <- nrow(segments)
	uni.segments <- as.data.frame(matrix(NA,ncol=7,nrow=0))
	colnames(uni.segments) <- c("sampleID",colnames(segments)[1:5],"mean")
	for(i in 1:nSample){
		sampleID <- colnames(segments)[5+i]
		sample.segments <- cbind(rep(sampleID,nSeg),segments[,c(1:5,(5+i))])
		colnames(sample.segments) <- colnames(uni.segments)
		uni.segments <- rbind(uni.segments,sample.segments,deparse.level=0)
	}

	return(uni.segments)
	
}