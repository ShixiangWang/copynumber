####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that makes it possible to plot frequencies based on segments data frame. 

### Required by:
### plotFreq
### plotCircle

### Requires:
### is.multiseg
### getUnisegFormat
### interpolate.pcf
### subsetSegments

getFreqData <- function(data){  
  #Check if segments or data:
  if(colnames(data)[1]=="sampleID" || colnames(data)[2]=="arm"){
    #input is segments data frame;
    #could be on a multiseg-format -> convert to uniseg-format:
    if(is.multiseg(data)){
      data <- getUnisegFormat(data)
    }
    #need to convert to appropriate format
    #first find intersection of all breakpoints:
    chr <- unique(data[,2])
    bpts <- matrix(NA,nrow=0,ncol=2)
    for(j in 1:length(chr)){
      subseg <- subsetSegments(data,chrom=chr[j])
      bpts <- rbind(bpts,data.frame(chr[j],sort(unique(c(subseg$start.pos,subseg$end.pos)))))
    } 
    colnames(bpts) <- c("chrom","pos")

    #Then interpolate to get pcf-val in all breakpoints
    data = interpolate.pcf(data, bpts) 
  }
  #If not segments, the input data should already be on an appropriate format (chrom, pos, estimates)
  
  return(data) 
} 