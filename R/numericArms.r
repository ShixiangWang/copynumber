####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that converts character arms to numeric


##Input:
### chrom : vector with chromosome numbers corresponding to each character arm
### char.arms : vector containing charcter arms; dentoed p or q

## Output:
### arms : vector with numeric arms calculated as chrom*2 - 1 or chrom*2

##Required by:
### adjustPos
### adjustSeg
### multipcf
### fastPcf
### pcf
### winsorize
### aspcf

##Requires:
###  none


numericArms <- function(chrom,char.arms){
  p.arm <- which(char.arms=="p")
  q.arm <- which(char.arms=="q")
  arms <- rep(NA,length(char.arms))
  arms[p.arm] <- chrom[p.arm]*2-1
  arms[q.arm] <- chrom[q.arm]*2
  
  return(arms)
}