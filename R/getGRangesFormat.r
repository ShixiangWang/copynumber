
getGRangesFormat <- function(segments){
  nc = ncol(segments)
  if(!is.multiseg(segments)){
    gr <- GRanges(seqnames=segments$chrom,
                ranges=IRanges(start=segments$start.pos,end=segments$end.pos,names=segments$sampleID))
    mcols(gr) <- DataFrame(segments[,c(3,6:nc)])            
  }else{
    gr <- GRanges(seqnames=segments$chrom,
                ranges=IRanges(start=segments$start.pos,end=segments$end.pos))
    mcols(gr) <- DataFrame(segments[,c(2,5:nc)])                   
  }
  
  return(gr)
  
}

