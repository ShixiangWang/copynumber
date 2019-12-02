
###################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# Function that separates chromosome arms in genome-plots by dashed lines

## Required by:
## plotFreq (genomeFreq)

## Requires:
## getArmandChromStop
## convert.unit



#Function used to separate chromosome arms by stapled lines in genome plot:
addArmlines <- function(chromosomes,xaxis,unit,ind=NULL,cex,op){
  if(xaxis=="pos"){
  	#Use cytoband data information to get stopping points of chromosome arms:
  	marks <- getArmandChromStop(op$assembly,unit)
  	armstop <- c(marks$pstop[1],cumsum(marks$chromstop)[1:length(marks$chromstop)-1]+marks$pstop[2:length(marks$pstop)])
  	scale.fac <- convert.unit(unit1=op$plot.unit,unit2=unit)    #Scaling factor according to plot.unit
    arm.mark <- armstop*scale.fac
    #Separate arms by vertical lines in existing plot:
    arg <- list(chrom.lwd=1, chrom.lty=2, chrom.col="darkgrey",chrom.side=3, chrom.cex=cex,chrom.line=c(0,0.3))
  
    if(!is.null(op)){
      arg <- modifyList(arg,op)
    }  
    abline(v=arm.mark[1:(length(arm.mark)-1)],col=arg$chrom.col,lwd=arg$chrom.lwd,lty=2)
  }
}


