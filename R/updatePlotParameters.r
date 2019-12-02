####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that updates pcfPlot parameters
#Updates or sets main, ylimits, xlimits, tickmarks, mgp, and margin

#Input:
### seg.lim: segmentation limits
### xmax: maximum on xaxis
### type: plot type
### xaxis: index or positions along xaxis?
### op: list with set plot parameters
### data.lim: data limits
### sampleID: the sampleID to be plotted
### k: the chromosome number to be plotted


###Output:
### op: list with updated plot parameters

### Required by:
### plotObs
### plotSegments

### Requires:
### get.xticks
### get.yticks


updatePlotParameters <- function(seg.lim,xmax,type,xaxis,op,data.lim,sampleID,k){

  #Set MAIN TITLE for this plot depending on type of plot:
  if(is.null(op$main)){
    if(type %in% c("genome","chromosome")){
      op$main <- sampleID
		}else{
		  #type = sample or aspcf
		  op$main <- paste("Chromosome",k,sep=" ")	
		}
	}#endif
	
  #YLIM;if not specified by user: 
  #leave out the q/2 % most extreme observations in both directions (either based on entire genome or within chromosome), make sure h is included in plot
  #and that segment-limits are taken into account:
  if(is.null(op$ylim)){
	   op$ylim[1] <- min(data.lim[1],op$h,seg.lim[1])  #y-min
	   op$ylim[2] <- max(data.lim[2],op$h,seg.lim[2])  #y-max
  }
  
	#set XLIM parameters:
  if(is.null(op$xlim)){
		op$xlim <- c(0,xmax)
	}


	#TICK MARKS on x and y axis:
	if(is.null(op$at.x)){
		if(type=="genome"){
			n.ticks <- 10
		}
		else{
			n.ticks <- 6
		}
		if(xaxis=="pos"){
			op$at.x <- get.xticks(op$xlim[1],op$xlim[2],unit=op$plot.unit,ideal.n=n.ticks)
		}else{
			#xaxis = index:
			op$at.x <- get.xticks(op$xlim[1],op$xlim[2],unit="mbp",ideal.n=n.ticks)
		}
	}
	if(is.null(op$at.y)){
		op$at.y <- get.yticks(op$ylim[1],op$ylim[2])
	}

  #MGP; placement of labels and axis annotation:
	if(is.null(op$mgp)){
		op$mgp <- c(1.3,0.05,0)*op$f
		mgp.y <- c(2,0.5,0)*op$f
	}else{
		mgp.y <- op$mgp
	}
  op$mgp.y <- mgp.y
  
  
  #Default is no xlab, if specified by user the BOTTOM MARGIN must be increased:
	if(op$xlab!=""){
		op$mar[1] <- op$mar[1] + 1
	}
  
  return(op)
}#end updatePlotParameters