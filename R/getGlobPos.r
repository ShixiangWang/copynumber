####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to covert local posistions to global positions based on a defined set of local stop positions chromosomes given in cyto.data

##Input:
### chromosomes: vector with chromosome numbers
### position: vector with positions
### pos.unit: positon's unit
### cyto.data: data frame with cytoband information
### delta: number of base pairs to add to end of each chromosome (used by plotCircle)

##Output:
### glob.pos: vector with global positions

##Required by:
### adjustSeg
### getx 
### plotCircle

##Requires:
### getArmandChromStop



getGlobPos <- function(chromosomes, position, pos.unit, cyto.data, delta=0){
  
	#Get local stopping posistions for each p-arm and each chromosome from cytoband data 
	l <- getArmandChromStop(cyto.data=cyto.data,unit=pos.unit)
  chromstop <- l$chromstop 
 
  #Need to make sure that none of the input positions are larger than chromosome stop postions:
  u.chrom <- unique(chromosomes)
  new.pos <- position
  for(j in 1:length(u.chrom)){
    probe.c <- chromosomes==u.chrom[j]
    #Check for positions that are larger than chromosome max position
    out <- which(position[probe.c] > chromstop[u.chrom[j]])
    #Replace these positions by max chrom position:
	  new.pos[probe.c][out] <- chromstop[u.chrom[j]] 
  }
  chromstop <- chromstop + delta
  glob.chromstop <- cumsum(chromstop)   #Global stopping position for each chromosome
  
  glob.pos <- new.pos + c(0,glob.chromstop[-length(glob.chromstop)])[chromosomes]  #Calculate global positions by adding global chromstop (for chrom > 1, for chrom=1 positions remain unchanged
	return(glob.pos)
	
}#end getGlobPos