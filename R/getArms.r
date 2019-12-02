####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Get arm numbers from chromosomes, position and cytoband data:

##Input:
### chrom: vector of chromosome numbers
### pos: vector of positions
### pos.unit: unit for positions
### cyto.data: object specifying which genome assembly version should be used for cytoband data

##Output:
### arms: a character vector with arms; dentoed p and q 

##Required by:
### adjustPos
### multipcf
### pcf
### winsorize
### aspcf 


##Requires:
### getArmandChromStop
### numericChrom

getArms <- function(chrom, pos, pos.unit="bp", cyto.data){

  #Make sure chromosomes are numeric:
  chrom <- numericChrom(chrom)
  
	nProbe <- length(chrom)
	chrom.list <- unique(chrom)
	nChrom <- length(chrom.list)
  
	#Get local stopping posistions for each p-arm and each chromosome from cytoband data 
	l <- getArmandChromStop(cyto.data=cyto.data,unit=pos.unit)
	pStop <- l$pstop
	chromStop <- l$chromstop 
	
	#Intitialize
	arms <- rep(NA,nProbe)
	
	for(i in 1:nChrom){
		#Find corresponding arm numbers:
		c <- chrom.list[i]
		ind.c <- which(chrom==c)
			
		arms[ind.c] <- "q"
		p.arm <- ind.c[pos[ind.c]<=pStop[c]]
		arms[p.arm] <- "p"   #p-arm
		
	}
	
	return(arms)
}

