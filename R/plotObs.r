####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

###Function that plots a vector of copy numbers

###INPUT:
### type: what kind of plot is this? may be "sample","chromosome","genome" or "aspcf"
### y   : a vector of data points to be plotted;  
### pos : a vector of positions to be plotted along x-axis
### unit: the unit used to represent probe positions (bp / kbp / mbp)
### winsoutliers: a vector giving winsorization outlier statuses for data points in y
### xaxis: a string indicating whether positions ("pos") or indeces ("index") should be plotted on xaxis
### plot.ideo: will ideogram be plotted?
### k: the chromosome number being plotted; only applicable for type="sample" or if print.xwarn..
### sampleID: the id of the sample being plotted; only applicable for type="chromosome" and type="genome"
### chromosomes: a vector of same length as pos giving the corrosponding chromosome numbers (only used when type="genome")
### frame: the frame dimensions for this plot
### new: has a new plot already been started?
### print.xwarn: should a warning be printed if the input positions exceed the chromosome range in cytoband?  st the moment always FALSE
### seg.lim:  ylimits for segments ; check whether this should be adjusted for q
### data.lim: ylimits for data points (adjusted for q and equalRange)
### op:  other parameters
### xmax: maximum on x-axis, will be NULL if ideogram is not plotted  


##Required by:
### plotAllele
### plotChrom
### plotGenome
### plotSample


##Requires:
### addChromlines
### addToPlot
### getPlotSymbols
### getx
### updatePlotParameters

plotObs <- function(type,y,pos,unit,winsoutliers,xaxis,plot.ideo=FALSE,k=NULL,sampleID=NULL,frame=NULL,new=FALSE,print.xwarn=FALSE,seg.lim=NULL,op,xmax=NULL,chromosomes=NULL,data.lim=NULL){


	#Pick out what should be plotted on x-axis (position or index). If type=genome positions are converted to global positions
  x <- getx(xaxis,type,chromosomes,pos,unit,op)
  if(is.null(xmax)){
    xmax <- max(x)
  } 

  #Update diverse plotparameters (ylim,xlim,main,ticks,mar)
  op <- updatePlotParameters(seg.lim=seg.lim,xmax=xmax,type=type,xaxis=xaxis,op=op,data.lim=data.lim,sampleID=sampleID,k=k)
  
  #Get colours, symbols and size for plotting observations. Truncate observations outside limits
  s <- getPlotSymbols(y,winsoutliers,type,x,xmax,print.xwarn,k,op)
  colobs <- s$colobs
  pch.obs <- s$pch.obs
  cex.obs <- s$cex.obs
  #The values for border-observations have been truncated
  x <- s$x
  y <- s$y   
	
	#PLOT OBSERVED DATA:

	#Empty plot with right dimensions:
	par(fig=unlist(frame),new=new,mar=op$mar)    
	plot(x,y,ylab="",xlab="",main="",pch=pch.obs,cex=cex.obs,col=colobs,ylim=op$ylim,xlim=op$xlim,xaxt="n",yaxt="n",xaxs="i",yaxs="r")

	#Add axes, labels, title and reference line
  addToPlot(plot.ideo,op,type=type)
	
	#Separate chromosomes by vertical lines (only done for type=="genome")
	if(type=="genome"){
    addChromlines(chromosomes,xaxis,unit,cex=op$cex.chrom,op=op)
	}
	


}#endplotObs




