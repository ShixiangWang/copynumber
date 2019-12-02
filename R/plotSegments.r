####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that plots segmentation results ; on top of existing plot or sets up a new plot

###INPUT:
### type: what kind of plot is this? may be "sample","chromosome","genome" or "aspcf"
### segments: a matrix/data frame with segments to be plotted;  
### unit: the unit used to represent probe positions (bp / kbp / mbp)
### xaxis: a string indicating whether positions ("pos") or indeces ("index") should be plotted on xaxis
### k: the chromosome number being plotted; seems only applicable for type="sample"
### sampleID: the id of the sample being plotted; only applicable for type="chromosome" and type="genome"
### frame: the frame dimensions for this plot
### plot.ideo: will ideogram be plotted?
### add: should segments plot be added on top of an existing plot?
### col,lty,lwd : color,line type and line width for this segmentation
### new: has a new plot been started?
### print.xwarn: should a warning be printed if the input positions exceed the chromosome range in cytoband? at the moment always FALSE
### seg.lim:  ylimits for segments ; check whether this should be adjusted for q
### data.lim: ylimits for data points
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
### adjustSeg
### connectSeg
### updatePlotParameters

plotSegments <- function(type,segments,unit,xaxis,k=NULL,sampleID=NULL,frame=NULL,plot.ideo=FALSE,add=FALSE,col,lty,lwd,new=FALSE,
                  print.xwarn=FALSE,seg.lim=NULL,data.lim=NULL,op,xmax=NULL){
	
	
  #Check that segments contains any information, and retrieve info
  if(nrow(segments)>0){
  
  	#Retrieve segmentinfo
  	nSeg <- nrow(segments)
  	chrom <- segments[,2]
  	arms <- segments[,3]
  	start <- segments[,4]
  	stop <- segments[,5]
  	nPos <- segments[,6]
  	seg.mean <- segments[,7]
  	
  	
  	#Adjust start and stop according to plot type, xaxis, and connect
  	a <- adjustSeg(chrom,arms,start,stop,nPos,type,xaxis,unit,op$connect,op)
    use.start <- a$use.start
    use.stop <- a$use.stop
    sep.arm <- a$sep.arm
    if(is.null(xmax)){	
      xmax <- max(use.stop)
    }
    
  
  	#Plotting:
  	if(add){
  		#add segments to existing plot:
  		segments(x0=use.start,y0=seg.mean,x1=use.stop,y1=seg.mean,col=col,lwd=lwd,lty=lty)
  		
  		#Should segments be connected? 
  		if(op$connect){
  		  connectSeg(sep.arm,nSeg,use.stop,seg.mean,col,lwd,lty)
  		}#endif
  		
  	}else{
      #Set up plot; only segments are to be plotted
     
      #Update diverse plotparameters (ylim,xlim,main,ticks,mar)
      op <- updatePlotParameters(seg.lim,xmax,type,xaxis,op,data.lim=data.lim,sampleID=sampleID,k=k)
        
  		#Check if end of any segments is larger than max in ideogram; if this is the case
  		#a warning is printed
  		out.seg <- use.stop[use.stop>xmax]		
  		if(length(out.seg)>0 && print.xwarn){
  			#Print warning:
  			warning(paste("Chromosome",k,"ranges from position 0 to",paste(xmax,".",sep=""),length(out.seg),"segments are outside this range.",sep=" "),
                    call.=FALSE,immediate.=TRUE)
  		}
  
  		#PLOTTING:
  		
  		#Empty plot with desired dimensions:
  		par(fig=unlist(frame),new=new,mar=op$mar)
      plot(1,1,type="n",main="",xlab="",
  			xlim=op$xlim,las=op$las,ylab="",ylim=op$ylim,xaxt="n",yaxt="n",xaxs="i")
  		 
  		#Plot segments
  		segments(x0=use.start,y0=seg.mean,x1=use.stop,y1=seg.mean,col=col,lwd=lwd,lty=lty)
  		
      #Should segments be connected? 
  		if(op$connect){
  		  connectSeg(sep.arm,nSeg,use.stop,seg.mean,col,lwd,lty)
  		}#endif
  		
  		#Add axes, labels, title and reference line
      addToPlot(plot.ideo,op,type=type)
		
  		#Separate chromosomes by vertical lines (only done for type=="genome")
      if(type=="genome"){
        addChromlines(chrom,xaxis,unit,ind=c(use.start,use.stop[length(use.stop)]),cex=op$cex.chrom,op=op)
      }
  		
  	}#endif
	
  }#endif	
}#endplotSegments




