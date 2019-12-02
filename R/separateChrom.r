####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that returns the index where each chromosome starts (and the last chromosome ends)

##Input:
### v: a vector of chromosome numbers

## Output:
### cp: indeces for start of each chromosome and end of last chromosome

##Required by:
### addChromlines
### adjustSeg

##Requires:
### none


separateChrom <- function(v){
	d <- diff(v)   #get difference between value (i+1) and value i in vector v
	cp <- which(d!=0)+1  #get changepoints
	
	#Add start of vector and (stop+1) of the whole vector
	cp <- c(1,cp,(length(v)+1))

	return(cp)
}#end separateChrom
