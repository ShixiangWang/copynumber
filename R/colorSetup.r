####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to set range of colors for heatplot:


##Input:
### upper.lim, lower.lim: limits for aberration calling
### op: list with other plot parameters

##Output:
### colors: a vector giving a range of color nuances
### intervals: vector giving intervals between lower.lim and upper.lim corresponding to different color nuances

##Required by:
### plotHeatmap 

##Requires:
### none

colorSetup <- function(upper.lim,lower.lim,op){
	#Colors:
	if(is.na(op$n.col)){
    op$n.col <- 50
	}
	#Range of colours in shades:
	range = upper.lim - lower.lim
  n.col1 = ceiling(op$n.col*abs(lower.lim)/range)
  n.col2 = ceiling(op$n.col*upper.lim/range)
  col.down = colorRampPalette(op$colors[1:2], space="rgb")(n.col1)
  col.up = colorRampPalette(op$colors[2:3], space="rgb")(n.col2)
  colors = c(col.down,col.up)
  
	#Divide into sequence according to values of limits:
	lim.seq1 <- seq(lower.lim,0,length.out=n.col1-1)
  lim.seq2 <- seq(0,upper.lim,length.out=n.col2-1)
  lim.seq <- c(-Inf,unique(c(lim.seq1,lim.seq2)),Inf)
  n.seq <- length(lim.seq)

	#Get matrix with intervals:
	low <- lim.seq[1:(n.seq-1)]
	high <- lim.seq[2:n.seq]
	intervals <- matrix(c(low,high),ncol=2,nrow=n.seq-1,byrow=FALSE)

	return(list(colors=colors,intervals=intervals))
}

