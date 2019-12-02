####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to set default plotting parameteres for frequency plots, and modify these according to user specifications:

#Input:
### type: plot type (genome or bychrom)
### nc: number of columns in plot
### nr: number of rows in plot
### thres.gain,thres.loss: aberration calling thresholds
### chrom: a vector giving the chromosomes to be plotted, only used for type=bychrom
### ... : other optional plot parameters specified by user

#Output:
### op: a list containing default and user modified plot parameters

###Required by: 
### plotFreq (genomeFreq and chromosomeFreq)
### plotWeightedFreq (weightedGenomeFreq and weightedChromosomeFreq)



getFreqPlotParameters <- function(type,nc,nr,thres.gain,thres.loss,chrom=NULL,...){

  #Apply a scaling factor according to number of columns and rows in plot:
	#seems to work ok:
	cr <- nc*nr
	f <- 1-0.013*cr


  #Common default parameters for genome and bychrom:
  op <- list(ylab="% with gain or loss",
             plot.size=c(11.8,min(3*nr,8.2)),
             col.gain="red",
             col.loss="blue",
             plot.unit="mbp",
             percentLines=TRUE,
             continuous=TRUE,
             assembly="hg19",
             las=1,
             chrom.lty=5,
             chrom.side=1,
             chrom.col="darkgrey",
             cyto.text=FALSE,
             #Parameters that will be set later
             ylim=NULL,
             xlim=NULL,
             xlab=NULL,
             mgp=NULL,
             mar=NULL,
             at.x=NULL,
             at.y=NULL,
             ideo.frac=NA,
             #Parameters that depend on the grid layout of the plot
             f=f,
             main.line=0.6*f,
             cex.lab=0.9*f,
             cex.main=f,
             cex.axis=0.8*f,
             cex.cytotext=0.6*f,
             cex.chrom=0.8*f
             )
  
  #Defaults specific to plot type:
  #For genome plot:
	if(type=="genome"){
		op$main <- paste("Thresholds = [",thres.loss,",",thres.gain,"]",sep="")
		op$main.line <- 1.5*f
    op$xlab <- ""
    op$plot.ideo <- FALSE
    op$mar <- c(1.5*op$f,3*op$f,2.3*op$f,1*op$f)
  }

	#For chromosome plot:
	if(type=="bychrom"){
		op$main <- paste("Chromosome ",chrom,sep="")
    op$title <- paste("Thresholds = [",thres.loss,",",thres.gain,"]",sep="")
    op$plot.ideo=TRUE
	}
	
	#Check for user modifications
	op <- modifyList(op,list(...))

  #Set/modify parameters more depending on user input:
  
  #Set assembly to refer to stored data instead of character string:
  op$assembly <- get(op$assembly)
  
  #Placement of labels and axis annotation:
	if(is.null(op$mgp)){
		op$mgp <- c(1.3,0.05,0)*f
		mgp.y <- c(2,0.5,0)*f
	}else{
		mgp.y <- op$mgp
	}
  op$mgp.y <- mgp.y

  #Xlabel
  if(is.null(op$xlab)){
    op$xlab <- paste("Position (",op$plot.unit,")",sep="")
  }

  #margins:
  if(is.null(op$mar)){
    op$mar <- if(op$plot.ideo) c(0.2*f,3*f,2.5*f,f) else c(1.5*f,3*f,2.5*f,f)    
  }
  
  #Set default ideo.frac and ideogram margins:
  if(op$plot.ideo){
    #ideogram margins:
    op$mar.i <- c(0.2*f,3*f,0,f)
    if(op$cyto.text){
      #Need to increase bottom margin:
      op$mar.i <- op$mar.i + c(2,0,0,0)
    }
    #Make sure left and right margins are equal for mar and mar.i:
    op$mar.i[c(2,4)] <- op$mar[c(2,4)] 
    if(is.na(op$ideo.frac)){
      #ideo.frac has not been defined by user:
      op$ideo.frac <- 0.05*sqrt(sqrt(cr))
      if(op$cyto.text){
        #Need larger space for ideogram:
        op$ideo.frac <- op$ideo.frac*2
      }
    }
  }else{
    op$ideo.frac <- 0
  } 

  #Check that we have a title for each plot:
  if(type=="genome" && length(op$main) < length(thres.gain)){
		op$main <- rep(op$main[1],length(thres.gain))
	}
	if(type=="bychrom" && length(op$main)<length(chrom)){
	 op$main <- rep(op$main[1],length(chrom))
	}
	#Make sure there is enough titles:
	if(type=="bychrom" && length(op$title) < length(thres.gain)){
		op$title <- rep(op$title[1],length(thres.gain))
	}
	
	return(op)
}#end function