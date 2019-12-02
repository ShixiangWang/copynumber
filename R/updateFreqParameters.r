####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################


#Function to set default ylim and at.y for plotFreq

##Input: 
### freq.del: vector with deletion frequencies
### freq.amp: vector with amplification frequencies
### op: list with plot parameters

##Output:
### op: list wiht updated plot parameters

##Required by: 
### plotFreq (genomeFreq and chromosomeFreq)
### plotWeightedFreq (weightedGenomeFreq and weightedChromosomeFreq)

##Requires: none


updateFreqParameters <- function(freq.del,freq.amp,op){
  #Y-limits; symmetric:
  max.freq <- max(c(freq.del,freq.amp))
  
  #Define tickmarks on y-axis
  if(max.freq>30){
    at.y <- seq(0,100,by=25)
  }else if(max.freq>10){
    at.y <- seq(0,100,by=10)
  }else{
    at.y <- seq(0,100,by=5)
  }
  if(is.null(op$at.y)){
    op$at.y <- c(at.y)
  }
  #Make sure ylim includes the first tickmark above max.freq
  q <- min(op$at.y[op$at.y>=max.freq])
  ylim <- c(-q,q)
    
  if(is.null(op$ylim)){
    op$ylim <- ylim
  }
  
  return(op)
}


  