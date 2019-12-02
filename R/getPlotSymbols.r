####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that sets colors, symbols and size for data points to be plotted
#Also truncates data to fit in plot

#Input:
### y: a vector with data points
### winsoutliers: a vector with outliers statuses for each observation in y
### type: plot type; sample, chromosome, genome or aspcf
### x: a vector with xaxis positions
### xmax: the maximum position
### print.warn: should a warning be printed if any positions are outside maximum chromosome positions?
### k: the chromosome being plotted
### op: other plot parameters

###Output:
### colobs: vector with colors to be used for each data point (winsorized data will have a different color)
### cex.obs: vector with plotsize for each data point
### pch.obs: vector with plotting symbol for each data point
### y: vector with data points; observations outside ylimits are truncated to fit in plot
### x: vector with x points; points outside xlimits are truncated to fit plot
 
##Required by:
### plotObs


##Requires:
###

getPlotSymbols <- function(y,winsoutliers,type,x,xmax,print.xwarn,k,op){
  #Colours, symbols and size for points to be plotted:
	n.k <- length(y)
	colobs <- rep(op$col,n.k)
	pch.obs <- rep(op$pch,n.k)
	cex.obs <- rep(op$cex,n.k)

	#Plot observations that fall outside y-range on the borders of the plot and mark by different symbol:
	small <- which(y < op$ylim[1])
	big <- which(y > op$ylim[2])
	y[small] <- op$ylim[1] 
	y[big] <- op$ylim[2] 
	
	colobs[c(small,big)] <- op$q.col
	cex.obs[c(small,big)] <- op$cex*op$q.cex
	pch.obs[c(small,big)] <- op$q.pch

  #Winsorized obs are plotted by different color:
	if(!is.null(winsoutliers)){
    outliers <- which(winsoutliers!=0)
		
		colobs[outliers] <- op$wins.col
		pch.obs[outliers] <- op$wins.pch
		cex.obs[outliers] <- op$cex*op$wins.cex
	}

	#Check if maximum probe position in data is larger than max in ideogram; if this is the case
	#a warning is printed
	out.pos <- x[x>xmax]
	if(length(out.pos)>0 && print.xwarn){
		#Print warning:
		warning(paste("Chromosome",k,"ranges from position 0 to",xmax,"mbp.",length(out.pos),"probe positions are outside this range.",sep=" "),call.=FALSE,immediate.=TRUE)
	}

	#Plot observations that fall outside x-range on the borders of the plot and mark by different symbol:
	cex.obs[x>xmax] <- op$cex*op$q.cex
	pch.obs[x>xmax] <- op$q.pch
	x[x>xmax] <- xmax

	return(list(colobs=colobs,cex.obs=cex.obs,pch.obs=pch.obs,y=y,x=x))
}#end getPlotSymbols