####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to set default plotting parameteres for plotsHeatmap, and modify these according to user specifications:

#Input:
### type: plot type (genome or bychrom)
### nc: number of columns in plot
### nr: number of rows in plot
### nSample: number of samples to include in plot
### upper.lim,lower.lim: aberration limits
### chrom: a vector giving the chromosomes to be plotted, only required if type=bychrom
### ... : other optional plot parameters specified by user

#Output:
### op: a list containing default and user modified plot parameters

###Required by: 
### plotHeatmap (genomeHeat and chromosomeHeat)


getHeatParameters <- function(type,nc,nr,nSample,upper.lim,lower.lim,chrom=NULL,ab=FALSE,...){

	#Apply a scaling factor according to number of columns and rows in plot:
	#virker som om dette funker ok:
	cr <- nc*nr
	f <- 1-0.013*cr

  
	
	#Common defaults for genome and chromosome heatplot:
  op <- list(colors=c("dodgerblue", "black", "red"),
             assembly="hg19",
             sep.samples=0,
             chrom.col="white",
             chrom.lty=2,
             chrom.side=1,
             sample.labels=TRUE,
             continuous=TRUE,
             plot.unit="mbp",
             sample.line=0.2,
             cyto.text=FALSE,
             cex.cytotext=0.6,
             #Parameters to be set later:
             ylab="",
             xlab=NULL,
             xlim=NULL,
             mar=NULL,
             plot.size=NA,
             n.col=NA,
             ideo.frac=NA,
             #Parameters that depend on the layout of the plot
             f=f,
             cex.main=f,
             cex.lab=f,
             cex.axis=0.8*f,
             main.line=0.6*f,
             sample.cex=0.7*f,
             cex.chrom=0.8*f
             )
   
  if(ab){
    op$colors <- c("dodgerblue","red")
    op$chrom.col <- "darkgrey"
    op$sep.samples <- 2/nSample
  }   
  #Additional defaults according to type
  if(type=="genome"){
    #Set default plot window size:
    op$plot.size[1] <- 11.8
    op$plot.size[2] <- ifelse(nr==1,min(nSample*0.8,8.2),min(8.2,3*nr))
    op$main <- paste("Limits = [",lower.lim,",",upper.lim,"]",sep="")
    op$main.line <- 1.7*f
    op$xlab <- ""
    op$mgp <- c(1.3,0.5,0)
    op$plot.ideo <- FALSE 
    op$mar <- c(1.5*op$f,3*op$f,2.5*op$f,1*op$f)
	}

	if(type=="bychrom"){
    #Set default plot size:
    op$plot.size[2] <- min(nSample*0.5*nr,8.2)
		op$plot.size[1] <- 11.8
    op$title <- paste("Limits = [",lower.lim,",",upper.lim,"]",sep="")
		op$main <- paste("Chromosome",chrom,sep=" ")
    op$mgp <- c(1.3,0.1,0)*f
    op$plot.ideo=TRUE
  }

  #Check for user modifications
	op <- modifyList(op,list(...))

  #Set assembly to refer to stored data instead of character string:
  op$assembly <- get(op$assembly)
  
  #More defaults/adjustments according to user input:
  
  #Xlabel
  if(is.null(op$xlab) && !op$plot.ideo){
    op$xlab <- paste("Position (",op$plot.unit,")",sep="")
  }
  
  #margins:
  if(is.null(op$mar)){
    op$mar <- if(op$plot.ideo) c(0.5*f,2*f,3*f,f) else c(2*f,2*f,3*f,f)
    if(op$sample.labels){
      op$mar[2] <- op$mar[2] + 2*op$f
    }
  }
  
  #Set default ideo.frac and ideogram margins:
  if(op$plot.ideo){
    #ideogram margins:
    op$mar.i <- c(0.2*f,2*f,0,f)
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
  if(type=="genome" && length(op$main) < length(upper.lim)){
		op$main <- rep(op$main[1],length(upper.lim))
	}
	if(type=="bychrom" && length(op$main)<length(chrom)){
	 op$main <- rep(op$main[1],length(chrom))
	}
	#Make sure there is enough titles:
	if(type=="bychrom" && length(op$title) < length(upper.lim)){
		op$title <- rep(op$title[1],length(upper.lim))
	}
  return(op)
  
}#end getHeatParameters
