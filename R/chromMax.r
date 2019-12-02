####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Get the maximum position on a given chromosome

#Input:
### chrom: number between 1-24 indicating which chromosome's ideogram is to be plotted
### cyto.data: data frame with cytoband-information
### unit: the unit used for positions in the plot
### cyto.unit: the unit used to represent positons in cyto.data

# Output:
### max.pos: a scalar giving the maximum position on this chromosome (given the cytoband information)

##Required by:
### plotFreq (chromosomeFreq)
### plotAllele
### plotChrom
### plotSample
### plotHeatmap
### plotWeightedFreq


## Requires:
### convert.unit


chromMax <- function(chrom,cyto.data,pos.unit,cyto.unit="bp"){
	
	#Get scaling factor for positions
	s <- convert.unit(unit1=pos.unit,unit2=cyto.unit)
	
  #Get the rows in cytoband data that correspond to this chromosome
	if(chrom==23){
    txt <- "chrX"
  }else if(chrom==24){
    txt <- "chrY"
  }else{
    txt <- paste("chr",chrom,sep="") 
  }
  rows <- which(cyto.data[,1]==txt)
  
  ##Get max position for this chromosome and scale according to pos.unit
  max.pos <- max(cyto.data[rows,3])
  max.pos <- max.pos*s
  
	return(max.pos)
	
}