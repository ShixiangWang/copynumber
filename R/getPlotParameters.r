####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to set default plotting parameteres for pcfPlots, and modify these according to user specifications:

#Input:
### type: plot type, genome, sample, chromosome or aspcf
### cr: number of columns times number of rows in grid layout of plot
### nSeg: number of segmentation results to be plotted
### sampleID: a vector giving the sampleID for each sample to be plotted; only applicable for type=sample and type=aspcf
### chrom: a vector giving the chromosomes to be plotted; only applicable for type=chromosome
### plot.ideo: is chromosome ideogram to be plotted
### xaxis: what is to be plotted along xaxis; may be "pos" or "index"
### ... : other optional plot parameters specified by user

#Output:
### arg: a list containing default and user modified plot parameters

##Required by:
### plotAllele
### plotChrom
### plotGenome
### plotSample




getPlotParameters <- function(type,cr,nSeg,sampleID=NULL,chrom=NULL,plot.ideo,xaxis,assembly,...){

	#Apply a scaling factor according to number of columns and rows in plot:
	#seems to work ok:
	f <- 1-0.013*cr
	
	
	#List with default plotting parameters:
	arg <- list(title="", 
              dir.print = NULL, 
              file.name = NULL, 
              onefile = TRUE,
              col="grey",
              wins.col="magenta",
              q.col="grey",
              pch=46,
              wins.pch=42,
              q.pch=42,
              ylab="Log R",
              las=1,
              h=0,
              h.lty=5,
              equalRange=TRUE,
              h.col="darkgrey",
              cyto.text=FALSE,
              plot.unit="mbp",
              seg.col=rainbow(nSeg),
              seg.lty=rep(1,nSeg),
              connect=TRUE,
              equalRange=TRUE,
              legend=ifelse(nSeg>1,TRUE,FALSE),
              plot.size=c(11.6,8.2),
              q=0.01,
              #Parameters that depend on the number grid layout in plot:
              f=f,
              mar=if(plot.ideo) c(0.2*f,3*f,2.5*f,f) else c(1.5*f,3*f,2.5*f,f),
              cex=2.5*f,
              wins.cex=0.4,
              q.cex=0.4,
              cex.lab=0.9*f,
              cex.main=f,
              cex.axis=0.8*f,
              cex.cytotext=0.7*f,
              cex.chrom=0.8*f,
              main.line=0.6*f,
              tcl=-0.3*f,
              h.lwd=2*f,
              seg.lwd=rep(3.5*f,nSeg),
              #Parameters that are set at later stage:
              xlab=NULL,
              main=NULL,
              xlim=NULL,
              ylim=NULL,
              at.x=NULL,
              at.y=NULL,
              mgp=NULL,
              ideo.frac=NA,
              assembly=assembly
              )


	if(type=="genome"){
		arg$main.line=1.7*f
		arg$xlab = ""
 }
	
 if(type=="sample"){
		arg$title <- sampleID
 }
 if(type=="chromosome"){
		arg$title <- paste("Chromosome",chrom,sep=" ")
 } 
	
	if(type=="aspcf"){
	 arg$ylab <- c("logR","BAF")
	 arg$h <- c(0,0.5)
	 arg$title <- sampleID
	}
	
	
	#Check for USER MODIFICATIONS:
	arg <- modifyList(arg,list(...))

  #Set assembly to refer to stored data instead of character string:
  arg$assembly <- get(arg$assembly)
  
  #Make sure ideogram is not plotted if xlim is specified:
  if(!is.null(arg$xlim)){
    plot.ideo = FALSE
  }
  
	#X-AXIS labels:
  if(is.null(arg$xlab)){
    if(xaxis=="index"){
      arg$xlab <- "Probe index"
    }else if(plot.ideo){
      arg$xlab <- ""
    }else{
      arg$xlab <- paste("Local position (",arg$plot.unit,")",sep="")
    }
  }#endif
	
	#Modify segment legend:
  if(is.logical(arg$legend)){
    if(!arg$legend){
      arg$legend <- NULL
    }else{
      arg$legend <- paste("Seg",1:nSeg,sep="")
    }
  }else{
    #Check length of user input for legend:
    if(!length(arg$legend)==nSeg){
      warning("Length of 'legend' does not match the number of segmentations, default legends are used instead",call.=FALSE)
      arg$legend <- paste("Seg",1:nSeg,sep="")
    }
  }#endif
  
  #Set default ideo.frac and ideogram margins:
  if(plot.ideo){
    #ideogram margins:
    arg$mar.i <- c(0.2*f,3*f,0,f)
    if(arg$cyto.text){
      #Need to increase bottom margin:
      arg$mar.i <- arg$mar.i + c(2,0,0,0)
    }
    #Make sure left and right margins are equal for mar and mar.i:
    arg$mar.i[c(2,4)] <- arg$mar[c(2,4)] 
    if(is.na(arg$ideo.frac)){
      #ideo.frac has not been defined by user:
      arg$ideo.frac <- 0.05*sqrt(sqrt(cr))
      if(arg$cyto.text){
        #Need larger space for ideogram:
        arg$ideo.frac <- arg$ideo.frac*2
      }
    }
  }else{
    arg$ideo.frac <- 0
  } 
  
  
	return(arg)
	
}#endgetPlotParameters
