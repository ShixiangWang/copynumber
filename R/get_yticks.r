
####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to calculate pretty ticks along y-axis in plot:

##Input:
### min: minimum observed y-value
### max: maximum observed y-value

##Output:
### ticks: tickmarks to be used along yaxis

##Required by:
### plotGamma
### updatePlotParameters

##Requires:
### none

get.yticks <- function(min,max){
	
	ideal.n <- 5
	by <- c(0.05,0.1,0.2,0.5,1,2,5,10,50,100,500,1000,2000,5000,10000)
	use.min <- rep(NA,length(by))
	use.max <- rep(NA,length(by))
	n.tick <- rep(NA,length(by))

	for(i in 1:length(by)){
		use.max[i] <- max
		if(max%%by[i]!=0){
			use.max[i] <- max + (by[i]-max%%by[i]) 
		}
		use.min[i] <- min
		if(min%%by[i]!=0){
			use.min[i] <- min - min%%by[i] 
		}
		seq <- seq(use.min[i],use.max[i],by=by[i])
		n.tick[i] <- length(seq)
				
	}#endfor	
		
	diff <- sapply(n.tick,"-",ideal.n)
	best <- which.min(abs(diff))
		
	ticks <- seq(use.min[best],use.max[best],by=by[best])
	return(ticks)

}#endfunction
