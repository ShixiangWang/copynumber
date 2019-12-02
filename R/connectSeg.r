####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Connect segments by vertical plotting. For sample and chromosome plots (sep.arm will be != NULL), segments are not connected across different arms.

##Input:
### sep.arm: vector giving index at which chromosome arms changes
### nSeg: number of segments
### use.stop: segment stopping positions
### seg.mean: mean value of segments
### col: segment color
### lwd: segment line width
### lty: segment line type


##Required by:
### plotGamma
### plotSegments

##Requires:
## none


connectSeg <- function(sep.arm,nSeg,use.stop,seg.mean,col,lwd,lty){
  if(!is.null(sep.arm)&&nSeg>2){
      #No connection across arms in the case of sample- or chromosomePlot, and across not adjecent chromosomes in the case of genomePlot
			 segments(x0=use.stop[-c((sep.arm-1),nSeg)],y0=seg.mean[-c((sep.arm-1),nSeg)],x1=use.stop[-c((sep.arm-1),nSeg)],y1=seg.mean[-c(1,sep.arm)],col=col,lwd=lwd,lty=lty)
  }else{if(is.null(sep.arm)&&nSeg>1){
			 segments(x0=use.stop[1:(nSeg-1)],y0=seg.mean[1:(nSeg-1)],x1=use.stop[1:(nSeg-1)],y1=seg.mean[2:nSeg],col=col,lwd=lwd,lty=lty)
  }}
}#end connectSeg