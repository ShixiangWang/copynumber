
####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

##Required by:
## pcf

## Requires:
## none 


#find nearest non-missing neighbour, and add missing obs to segment where nearest neighbour is located
handleMissing <- function(nn,pos,obs,pos.start,pos.stop,seg.npos){

  #find out which segment interval nn belongs to, 
  inInt <- findInterval(pos[nn],pos.start)
  
  #add missing pos to number of pos in this segment 
  seg.npos[unique(inInt)] <- as.vector(seg.npos[unique(inInt)] + table(inInt))
  
  #move segment start or stop if missing obs pos is smaller or larger than existing segment
  pos.na <- pos[!obs]
  t1 <- which(pos.na < pos.start[inInt])
  pos.start[inInt[t1]] <- pos.na[t1]
  t2 <- which(pos.na > pos.stop[inInt])
  pos.stop[inInt[t2]] <- pos.na[t2]
  
  return(list(pos.start=pos.start,pos.stop=pos.stop,seg.npos=seg.npos))
  
}