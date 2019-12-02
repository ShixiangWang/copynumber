####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#function that find global limits on xaxis (used in genome plots)

##Input: 
### op: list with plot parameters
### pos.unit: the unit used for positions in data
### chrom: a vector of unique chromosome numbers found in data

##Output:
### op: a list with updated plot parameters (xlim)

##Required by:
### plotFreq (genomeFreq)
### plotHeatmap (genomeHeat)
### plotWeightedFreq
### plotGenome


##Requires:
### getArmandChromStop
### convert.unit


getGlobal.xlim <- function(op,pos.unit,chrom){  
  #Set xlim using chromosome information in cytoband; must transform to global information
  chromstop <- getArmandChromStop(op$assembly,pos.unit)$chromstop
  glob.chromstop <- cumsum(chromstop)   #Global stopping position for each chromosome
  scale.fac <- convert.unit(unit1=op$plot.unit,unit2=pos.unit)    #Scaling factor according to plot.unit
  #Not use chromosome X and Y if not in data
  if(!any(chrom==24)){
    glob.chromstop <- glob.chromstop[-24]
  }
  if(!any(chrom==23)){
    glob.chromstop <- glob.chromstop[-23]
  }
  xlim <- c(0,glob.chromstop[length(glob.chromstop)])*scale.fac
  
  return(xlim) 
}  