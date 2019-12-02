####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################


# Function to calculate running median for a given a window size

##Input:
### x: vector of numeric values
### k: window size to be used for the sliding window (actually half-window size)

## Output:
### runMedian : the running median corresponding to each observation

##Required by:
### getMad
### medianFilter


##Requires:
### none

medianFilter <- function(x,k){
  n <- length(x)
  filtWidth <- 2*k + 1
  
  #Make sure filtWidth does not exceed n
  if(filtWidth > n){
    if(n==0){
      filtWidth <- 1
    }else if(n%%2 == 0){
      #runmed requires filtWidth to be odd, ensure this:
      filtWidth <- n - 1
    }else{
      filtWidth <- n
    }
  }
  
  runMedian <- runmed(x,k=filtWidth,endrule="median")

  return(runMedian)

}
