####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

##Input:
### chromosomes: a vector of chromosome numbers
### xaxis: what should be plotted along xaxis; pos or index
### unit: the unit used for positions, bp, kbp or mbp
### ind: only used when called by plotSegments; vector with start position and last stop position for segments
### cex: size used to plot chromosome numbers
### op: other plot parameters

##Output:
### adds chromosome lines to existing genome plot

##Required by:
### plotObs
### plotSegments
### plotFreq
### plotHeatmap
### plotWeightedFreq


##Requires:
### convert.unit
### getArmandChromStop
### separateChrom

#Function used to separate chromosomes by stapled lines in genome plot:
addChromlines <- function(chromosomes,xaxis,unit,ind=NULL,cex,op){
  if(xaxis=="pos"){
			#Use cytoband data information to get stopping points of chromosomes:
			chromstop <- getArmandChromStop(op$assembly,unit)$chromstop
			scale.fac <- convert.unit(unit1=op$plot.unit,unit2=unit)    #Scaling factor according to plot.unit
			chrom.mark <- c(1,cumsum(chromstop))*scale.fac
			#Drop chrom24 if no observations for this chrom in data/segments:
			if(any(chromosomes==24)){
			 chromosomes <- 1:24
			}else{
			 chromosomes <- 1:23
			 chrom.mark <- chrom.mark[-length(chrom.mark)]
			}
  }else{
		  chrom.mark <- separateChrom(chromosomes)
		  if(!is.null(ind)){
		    #Segments are used to decide where to separate chromosomes; get index start where chromosome number changes and last index
		    chrom.mark <- ind[chrom.mark]
		  }
		  chrom.mark <- chrom.mark - 0.5
  }
	#Separate chromosomes by vertical lines in existing plot:
	nChrom <- length(unique(chromosomes))
	arg <- list(chrom.lwd=1, chrom.lty=2, chrom.col="darkgrey",chrom.side=3, chrom.cex=cex,chrom.line=c(0,0.3))
	
	if(!is.null(op)){
		arg <- modifyList(arg,op)
	}
	abline(v=chrom.mark[2:(length(chrom.mark)-1)],col=arg$chrom.col,lwd=arg$chrom.lwd,lty=arg$chrom.lty)
	
	at <- (chrom.mark[1:nChrom]-1)+(chrom.mark[2:(nChrom+1)]-chrom.mark[1:nChrom])/2
	chrom.names <- unique(chromosomes)
	#Plot half at bottom, half at top:
	bot <- seq(1,length(chrom.mark),2)
	top <- seq(2,length(chrom.mark),2)
  mtext(chrom.names[bot],side=1,line=arg$chrom.line[1],at=at[bot],cex=arg$chrom.cex)
	mtext(chrom.names[top],side=3,line=arg$chrom.line[2],at=at[top],cex=arg$chrom.cex)
	
}#endaddChromlines
