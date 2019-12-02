
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
### nPos:  vector giving the number of probes in the segments
### type:  plot type
### xaxis: pos or index
### unit: unit used to represent positions
### connect: should segments be connected (meet halfway between probes? )
### op: list with other plot parameters

##Output:
### use.start: adjusted segment start to be used in plot
### use.stop:  adjusted segement stop to be used in plot
### sep.arm: vector giving the index at which arm numbers are changed

##Required by:
### plotGamma
### plotSegments


##Requires:
### convert.unit
### getGlobPos
### numericArms
### separateChrom
### getArmandChromStop
### numericChrom


adjustSeg <- function(chrom,char.arms,start,stop,nPos,type,xaxis,unit,connect,op){

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

  #What  should be plotted along x-axis
	if(xaxis=="pos"){
		#Plot against probe position on x-axis

		if(type=="genome"){
			#Use global positions
      start <- getGlobPos(chrom,start,pos.unit=unit,cyto.data=op$assembly)
      stop <- getGlobPos(chrom,stop,pos.unit=unit,cyto.data=op$assembly)
      #Retrieve chromosomestop positions from cytoband data
      chromstop <- getArmandChromStop(cyto.data=op$assembly,unit=unit)$chromstop
      glob.chromstop <- cumsum(chromstop)
    }

		use.start <- start
		use.stop <- stop

    #Connect segments halfway between end and start of two segments
    if(connect){
		  if(nSeg>1){
			 #The segments should start and end halfway between two probes:
			 half <- (use.start[2:nSeg]-use.stop[1:(nSeg-1)])/2
			 use.start[2:nSeg] <- use.start[2:nSeg] - half
			 use.stop[1:(nSeg-1)] <- use.stop[1:(nSeg-1)]+half
			 
       #If type=genome: want connected segments for different chromosomes to occur at chromosome-lines; first segment in each chromosome will therefore start at the global chromosome start position and end at the global chromosome end position.
       if(type=="genome"){
        unik.chrom <- unique(chrom)
        use.stop[c(sep-1,nSeg)] <- glob.chromstop[unik.chrom]
        #use.start[c(1,sep)] <- c(0,glob.chromstop[unik.chrom[-length(unik.chrom)]]+1)
        use.start[c(1,sep)] <- c(use.start[1],glob.chromstop[unik.chrom[-1]-1]+1)

       }else{
        #No connection of segments across arms when type is "chromosome" or "sample"
        use.start[sep] <- start[sep]
        use.stop[sep-1] <- stop[sep-1]
       }

		  }#endif
    }
    #Want to scale the x-axis to fit the desired unit given in plot.unit (default is mega base pairs)
    scale.fac <- convert.unit(unit1=op$plot.unit,unit2=unit)
		use.start <- use.start*scale.fac
		use.stop <- use.stop*scale.fac
  }else{
		#xaxis=="index"
		#Plot against probe index on x-axis:

    use.start <- rep(NA,nSeg)
    use.start[1] <- 1
    use.stop <- rep(NA,nSeg)
    use.stop[nSeg] <- sum(nPos)

    if(nSeg>1){
      for(i in 2:nSeg){
				use.start[i] <- use.start[i-1] + nPos[i-1]
      }
      use.stop[1:(nSeg-1)] <- use.start[2:nSeg]-1

      if(connect){
        #The segments should start and end halfway between two probes:
        use.start[2:nSeg] <- use.start[2:nSeg]-0.5
        use.stop[1:(nSeg-1)] <- use.stop[1:(nSeg-1)]+0.5

        if(type!="genome"){
          #Start of first segment in second arm should not be halfway between probes
		      use.start[sep] <- use.start[sep]+0.5
		      use.stop[sep-1] <- use.stop[sep-1]-0.5
        }#endif
      }
    }
	
	}#endif
	
	#sep.arm is later used to connect segments. When type=genome we want to connect across chromosomes except when chromosome represented in segments are not 
  #adjecent. When type!=genome we do not want to connect across arms
  sep.arm <- NULL
  if(type=="genome"){
    unik.chrom <- unique(chrom)
    dchr <- diff(unik.chrom)
    sep.use <- which(dchr!=1)
    if(length(sep.use)!=0){
      sep.arm <- sep[sep.use]
    }
  }else{
    if(length(sep)!=0){
      sep.arm=sep
    }
  }

	return(list(use.start=use.start,use.stop=use.stop,sep.arm=sep.arm))
  
}#endfunction
