####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to find frame dimensions when a window is to be sectioned into nrow rows and ncol columns:

##Input:
## nrow: number of rows in plot
## ncol: number of columns in plot

##Output:
# a list giving the left, right, bottom and top dimensions for each of the nrow*ncol frames to be made in plot

##Required by:
### plotAllele
### plotChrom
### plotGenome
### plotSample
### plotFreq
### plotHeat
### plotWeightedFreq 


##Requires:
### none

framedim <- function(nrow,ncol){
	cl <- 0:(ncol-1)
	cr <- 1:ncol
	left <- rep(1/ncol,ncol)*cl
	right <- rep(1/ncol,ncol)*cr

	rt <- nrow:1
	rb <- (nrow-1):0
	top <- rep(1/nrow,nrow)*rt
	bot <- rep(1/nrow,nrow)*rb

	return(list(left=left,right=right,bot=bot,top=top))
}#endframedim
