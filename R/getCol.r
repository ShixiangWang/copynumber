####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to retrieve color in heatplot to be used for each segment:


##Input:
### seg.mean: the segment mean value
### colors: the range of color nuances to be applied
### intervals: the intervals of values corresponding to each color nuance

##Output:
### col: the color to be plotted for this segment

##Required by:
### plotHeatmap 

##Requires:
### none

#Get segments colors:
getCol <- function(seg.mean,colors,intervals){
	#Find the interval in which seg.mean is located:
	inInt <- which(intervals[,1]<seg.mean & intervals[,2]>=seg.mean)
	#Pick colour:
	col <- colors[inInt]
	return(col)
}