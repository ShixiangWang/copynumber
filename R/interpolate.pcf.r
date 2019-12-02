####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

### Requires:
### pullOutContent

#Function with interpolates pcf-segments
interpolate.pcf = function(segments, x) {
  
  #Make sure segments is a data frame
  segments <- pullOutContent(res=segments,what="segments")
  
  usamp = unique(segments$sampleID)
  nsamp = length(usamp) 
  chrom = unique(x[,1])
  z = data.frame(cbind(x[,c(1:2)],matrix(0,nrow(x),nsamp)))
  #z = data.frame(x[,c(1,2)], matrix(0, nrow(x), nsamp))
  names(z) = c("chr","pos",usamp)
  for (i in 1:nsamp) {
    for (j in 1:length(chrom)) {
      fitij = segments[segments$sampleID==usamp[i] & segments$chrom==chrom[j],]
      v = (c(fitij$start.pos[-1],10^9)+fitij$end.pos)/2
      xj = x[x[,1]==chrom[j],2]
      kj = rep(0,length(xj))
      for (k in rev(1:length(v))) {
        kj[xj <= v[k]] = k
      }
      z[z$chr==chrom[j],2+i] = fitij$mean[kj]
    }
  }
  return(z)
}


