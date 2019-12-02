####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that checks if chrom is numeric, converts x/X and y/Y to 23 and 24 if not:

##Input:
### chrom: vector with chromosomes; numeric or character

## Output:
### chrom: numeric vector with chromosome numbers

##Required by:
### getArms
### checkSegments
### multipcf
### selectSegments
### selectData
### fastPcf
### pcf
### checkAndRetrievePlotInput
### plotFreq
### plotHeatmap
### plotWeightedFreq
### winsorize
### aspcf
### adjustSeg
### checkChrom

##Requires:
### none


numericChrom <- function(chrom){ 
 if(!is.numeric(chrom)){
    if(is.factor(chrom)){
      #If chrom is factor; need to convert to character first
      chrom <- as.character(chrom)
    }
    #Replace X by 23:
    chrx <- c(which(chrom=="x"),which(chrom=="X"))
    chrom[chrx] <- 23
    #Replace Y by 24
    chry <- c(which(chrom=="y"),which(chrom=="Y"))
    chrom[chry] <- 24
    
    chrom <- as.numeric(chrom)
  }
  return(chrom)
}