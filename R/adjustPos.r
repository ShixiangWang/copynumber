####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that scales positions according to plotunit, converts to global positions if type=genome, and finds start and stop positions (left and right) for recatangles to be plotted. Also makes sure that frequencies are shown as continuous if this is desired

##Input:
### position: the genomic postions to be plotted 
### chromosomes: the chromosomes corresponding to the positions
### pos.unit: the unit used for positions (bp,kbp,mbp)
### type: plot type (genome or chromosome)
###op: a list of other set plot parameters 

##Output:
### xleft: the left/start position of the plot rectangle
### xright: the right/stop position of the plot rectangle

##Required by:
### plotFreq (genomeFreq and chromosomeFreq)
### plotWeightedFreq (weightedGenomeFreq and weightedChromosomeFreq)


##Requires:
### getx
### getArms
### numericArms


adjustPos <- function(position,chromosomes,pos.unit,type,op){

  if(type=="chromosome"){
    #Only need to scale positions first
    pos <- getx(xaxis="pos",type=type,chromosomes=NULL,pos=position,unit=pos.unit,op=op)
  }else if(type=="genome"){
    #Need to scale and convert to global pos:
    pos <- getx(xaxis="pos",type=type,chromosomes=chromosomes,pos=position,unit=pos.unit,op=op)
  }
  
  nPos <- length(position)
  #Define left-pos and right-pos for freqency-rectangle to be plotted:
  xleft <- pos
  xright <- pos
  
  #Should frequencies be plotted continously across probes?:
  if(op$continuous){
    #The rectangles should start and end halfway between two probes, except across different arms/chromosomes, fixing this below
    half <- (pos[2:nPos]-pos[1:(nPos-1)])/2
    xleft[2:nPos] <- xleft[2:nPos] - half
    xright[1:(nPos-1)] <- xright[1:(nPos-1)] + half   
  }else{
    #Let rectangle be one probe wide:
    xleft[2:nPos] <- xleft[2:nPos] - 0.5
    xright[1:(nPos-1)] <- xright[1:(nPos-1)] + 0.5
  }
  
  
  if(type!="genome"){
    #First find locations for change in arm number:
    char.arms <- getArms(chromosomes,position,pos.unit,op$assembly)
    arms <- numericArms(chromosomes,char.arms)
    #Where do arms start:
    n.arm <- length(unique(arms))
  
    if(n.arm>1){
      #Locate where arm starts (if more than one arm):
      sep.arm <- separateChrom(arms)
      sep.arm <- sep.arm[-c(1,length(sep.arm))]
      #Keep positions at arm-change 
      xleft[sep.arm] <- pos[sep.arm]
      xright[sep.arm-1] <- pos[sep.arm-1]
    }
  }else{
    #when type=genome: only want to prevent continous over chromosomes:
    n.chrom <- length(unique(chromosomes))
    if(n.chrom > 1){
      sep.chrom <- separateChrom(chromosomes)
      sep.chrom <- sep.chrom[-c(1,length(sep.chrom))]
      #keep positions at chrom change
      xleft[sep.chrom] <- pos[sep.chrom]
      xright[sep.chrom-1] <- pos[sep.chrom-1]   
    }
  }
  
  return(list(xleft=xleft,xright=xright))
  
}#end adjustFreqPlotPos