####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Returns the values to be plotted along xaxis depending on the choice of xaxis (index or pos), and the type of plot
#If type=genome and xaxis=pos, global positions are calculated.
#Positions are scaled according to plotunit

### input:
### xaxis: index or pos
### type: plot type; genome, chromosome, sample or aspcf
### chromosomes: vector of same length as pos giving the corresponding chromosome numbers
### pos: vector giving probe positions 
### unit: unit used to represent positons   (bp,kbp, or mbp)
### op: a list containing other plot parameters

###Output:
### x: a vector containg the numbers to be plotted along xaxis of plot


##Required by:
### adjustPos
### plotObs

##Requires:
### getGlobPos
### convert.unit


getx <- function(xaxis,type,chromosomes,pos,unit,op){
  if(xaxis=="pos"){

    x <- pos
		if(type=="genome"){
			#Convert to global position:
			global.pos <- getGlobPos(chromosomes,pos,pos.unit=unit,cyto.data=op$assembly)   
			x <- global.pos
		}
		
		#Want to scale the x-axis to fit the desired unit given in plot.unit (default is mega base pairs)
    scale.fac <- convert.unit(unit1=op$plot.unit,unit2=unit)
		x <- x*scale.fac

	}else{
		#xaxis=="index"
		x <- 1:length(pos) 
	}#endif

	return(x)
	
}#endgetx