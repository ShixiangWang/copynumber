
####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

# Function that adjusts segment start and stop positions according to plot type, plot unit, xaxis etc

##Input:
### chrom: vector with segments chromosome numbers
### char.arms: vector with character arms
### start: start positions of segments
### stop: stop positions of segments
### type:  plot type
### unit: unit used to represent positions
### op: list with other plot parameters

##Output:
### use.start: adjusted segment start to be used in plot
### use.stop:  adjusted segement stop to be used in plot

##Required by:
### plotHeat


##Requires:
### convert.unit
### getGlobPos
### numericArms
### separateChrom
### numericChrom


adjustSegPos <- function(chrom,char.arms,start,stop,type,unit,op){

  #Make sure chromosome is numeric:
  chrom <- numericChrom(chrom)
  nSeg <- length(chrom)
  
  #Get indeces where chromosome number or arm changes
  if(type=="genome"){
    sep <- separateChrom(chrom)
  }else{
    #Convert character arms to numeric:
    arms <- numericArms(chrom,char.arms)
    sep <- separateChrom(arms)
  }
  #separateChrom adds indeces 1 and nSeg; don't need these here:
  sep <- sep[-c(1,length(sep))]

	if(type=="genome"){
		#Use global positions
    start <- getGlobPos(chrom,start,pos.unit=unit,cyto.data=op$assembly)
    stop <- getGlobPos(chrom,stop,pos.unit=unit,cyto.data=op$assembly)
  }

	use.start <- start
	use.stop <- stop

  #Continuous heatmap? connect halfway between end and start of two segments
  if(op$continuous){
	  if(nSeg>1){
		 #The heatmap should start and end halfway between two probes, except across arms/chromosomes 
		 half <- (use.start[2:nSeg]-use.stop[1:(nSeg-1)])/2
		 use.start[2:nSeg] <- use.start[2:nSeg] - half
		 use.stop[1:(nSeg-1)] <- use.stop[1:(nSeg-1)]+half
		 
     use.start[sep] <- start[sep]
     use.stop[sep-1] <- stop[sep-1]

	  }#endif
  }
  #Want to scale the x-axis to fit the desired unit given in plot.unit (default is mega base pairs)
  scale.fac <- convert.unit(unit1=op$plot.unit,unit2=unit)
	use.start <- use.start*scale.fac
	use.stop <- use.stop*scale.fac	
	
	return(list(use.start=use.start,use.stop=use.stop))
}#endfunction
