####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that finds nearest non-missing neighbour

##Input:
#pos: a vector of positions 
#obs: a logical vector indication which of the positions have observed copy number values

##Output:
#NNM: a vector which gives the nearest non-missing neighbour of the missing observations

##Required by:
### fastPcf
### pcf
### aspcf


##Requires:
### none


findNN <- function(pos,obs){
	
	ind.obs <- which(obs)
	pos.obs <- pos[obs]
	pos.na <- pos[!obs]
	
	#Find distances between position of missing obs and positions of non-missing obs
	d <- sapply(pos.na,FUN="-",y=pos.obs)
	if(!is.matrix(d)){
    d <- matrix(d,nrow=length(pos.obs),ncol=length(pos.na))
	}
	d <- abs(d)
	nn <- apply(d,2,which.min)  #find which observed obs is the nearest
	nn <- ind.obs[nn]  #get index for the nearest

	return(nn)

}
