
####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Input:
### op: a list with plot parameters
### type: plot type (genome or bychrom)

##Required by:
### plotFreq (genomeFreq and chromosomeFreq)
### plotWeightedFreq (weightedGenomeFreq and weightedChromosomeFreq)

##Requires:
### get.xticks


##Function that adds percentagelines, yaxis, xaxis and labels to frequency plots 
addToFreqPlot <- function(op,type){
    
  #Add y-lines if wanted
  if(is.logical(op$percentLines)){
    if(!op$percentLines){
      op$percentLines<- NULL
    }else{
      op$percentLines <- op$at.y
    }
  }
  if(!is.null(op$percentLines)){
    abline(h=c(-op$percentLines,op$percentLines),lty=3,col="grey82")
  }
  
  #Add yaxis and lab:
  #Make sure tickmarks are at percentLines
  if(!is.null(op$percentLines)){
    op$at.y <- op$percentLines
  }
  op$at.y <- c(-op$at.y,op$at.y)
  axis(side=2,cex.axis=op$cex.axis,at=op$at.y,mgp=op$mgp.y,las=op$las,tcl=-0.2,labels=abs(op$at.y))
  title(ylab=op$ylab,cex.lab=op$cex.lab,line=op$mgp.y[1])
  
  #Add xaxis:
  if(op$plot.ideo || type=="genome"){
    axis(side=1,labels=FALSE,tcl=0)
  }else{
    if(is.null(op$at.x)){
		  op$at.x <- get.xticks(0,op$xlim[2],unit=op$plot.unit,ideal.n=6)
		}
    axis(side=1,tcl=-0.2,at=op$at.x,cex.axis=op$cex.axis,mgp=op$mgp)
		title(xlab=op$xlab,cex.lab=op$cex.lab,line=op$mgp[1])
  }
  
  
}#end addToFreqPlot