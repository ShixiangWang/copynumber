
####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

##Input:
### plot.ideo: should ideogram be plotted, TRUE/FALSE
### op: other plot parameters
### type: plot type; genome, sample, chromosome or aspcf



##Required by:
### plotObs
### plotSegments

##Requires:
###
###

#Add axes, labels, title and reference line to plot
addToPlot <- function(plot.ideo,op,type){
  if(!plot.ideo && type!="genome"){
	  #Labels along xaxis
		axis(side=1,cex.axis=op$cex.axis,at=op$at.x,labels=as.integer(op$at.x),mgp=op$mgp,tcl=op$tcl)
	}else{
	  #No labels
		axis(side=1,labels=FALSE,tcl=0,at=op$xlim)
	}
	axis(side=2,cex.axis=op$cex.axis,at=op$at.y,mgp=op$mgp.y,las=op$las,tcl=op$tcl)

	#x- and y-lab
	mtext(text=op$xlab,side=1,line=op$mgp[1],cex=op$cex.lab)
	mtext(text=op$ylab,side=2,line=op$mgp.y[1],cex=op$cex.lab)

	#main title for this plot
	title(main=op$main,line=op$main.line,cex.main=op$cex.main)

  #Add reference line at y=h; h=NULL suppresses plotting of ref.line:
	if(!is.null(op$h)){
		abline(h=op$h,lty=op$h.lty,lwd=op$h.lwd,col=op$h.col)
	}#endif

}#end addToPlot