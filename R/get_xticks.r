####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to generate nice ticks along x-axis in plot:

##Input:
### min: min observed x-value
### max: max observed x-value
### unit: the plot unit for positions, either mbp or kbp
### ideal.n: the ideal number of tickmarks

##Output:
### ticks: a vector of nice tick marks

##Required by:
### plotGamma
### updatePlotParameters
### addToFreqPlot
### plotHeatmap


##Requires:
### none

get.xticks <- function(min,max,unit,ideal.n){
	
	if(!unit%in%c("mbp","kbp")){
	 stop("plot.unit must be one of 'mbp' and 'kbp'",call.=FALSE)
	}
	if(identical(unit,"mbp")){
		by <- c(1,2,5,10,20,40,60,80,100,200,500,1000,2000,5000)
	}else{if(identical(unit,"kbp")){
		ideal.n <- ideal.n-1
		by <- c(1,2,5,10,20,40,60,80,100,200,500,1000,2000,5000)*1000
	}}
	
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
	#Eventuelt den som minimerer positiv diff?	
	ticks <- seq(use.min[best],use.max[best],by=by[best])
	
	return(ticks)

}#endfunction