####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to get a scaling factor such that positions may be converted from unit2 to unit1:

##Input:
### unit1: unit that we want to convert to
### unit2: unit that should be converted

##Output:
### factor: a scaling factor used in conversion of positions

##Required by:
### addChromlines
### adjustSeg
### plotGamma
### chromMax
### getx
### plotIdeogram
### getGlobal.xlim

##Requires:
#none 

convert.unit <- function(unit1,unit2){
	
	factor <- NA
	#Allowed units:
	units <- c("bp","kbp","mbp")
	
	if(identical(unit1,unit2)){
		factor <- 1
	}else{
		if(identical(unit1,units[3])){
			if(identical(unit2,units[2])){
				factor <- 1/(10^3)
			}else{if(identical(unit2,units[1])){
				factor <- 1/(10^6)
			}}
		}else{
			if(identical(unit1,units[2])){
				if(identical(unit2,units[3])){
					factor <- 10^3
				}else{if(identical(unit2,units[1])){
					factor <- 1/(10^3)
				}}
			}else{
				if(identical(unit1,units[1])){
					if(identical(unit2,units[3])){
						factor <- 10^6
					}else{if(identical(unit2,units[2])){
						factor <- 10^3
					}}
				}
			}
		}
	}
	
	if(is.na(factor)){
		if(all(units!=unit1)){
			stop("plot.unit must be one of 'kbp' and 'mbp'",call.=FALSE)
		}
		
	}
	
	return(factor)
}#end convert.unit

