####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# Function that separates chromosomes in genome-plots by color

## Required by:
## plotFreq (genomeFreq)

## Requires:
## getArmandChromStop
## convert.unit



chromPattern <- function(pos.unit,op) {
  #Use cytoband data information to get stopping points of chromosomes:
	chromstop <- getArmandChromStop(op$assembly,pos.unit)$chromstop
	scale.fac <- convert.unit(unit1=op$plot.unit,unit2=pos.unit)    #Scaling factor according to plot.unit
	chrom.mark <- c(1,cumsum(chromstop))*scale.fac
	#Drop chrom24 if no observations for this chrom in data/segments:
	#if(!any(segments[,2]==24)){
	# chrom.mark <- chrom.mark[-length(chrom.mark)]
	#} 
	#Let background be black to avoid white parts in arms without probes:
  #rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "white")
  for (i in 1:(length(chrom.mark)-1)) {
    if(i%%2==0){
       rect(chrom.mark[i], par("usr")[3], chrom.mark[i+1], par("usr")[4], col = "grey95")#, border=NA)
    } 
  }
}
