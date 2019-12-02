####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that checks and returns valid filename(s):

##Input:
### nPage: number of plot pages
### file.name: input desired file.name
### ID: a vector of IDs for each plot page; eg chromosome number or sample name
### type: plot type

##Output:
### file.name: a vector giving valid filename(s)

##Required by:
### plotAllele
### plotChrom
### plotGenome
### plotSample




##Requires:
### none

getFilename <- function(nPage,file.name,ID,type){
  
  def.name <- switch(type,
                    genome="genomeplot",
                    sample="sampleplot",
                    chromosome="chromosomeplot",
                    aspcf="aspcfplot",
                    genomeFreq="genome frequencyplot",
                    chromFreq="chromosome frequencyplot",
                    genomeHeat="genome heatplot",
                    chromHeat="chromosome heatplot")
  if(nPage==1){
    #Only one file/window
    if(is.null(file.name)){
      file.name <- ifelse(length(ID)==1,paste(def.name,"_",ID,sep=""),def.name)
    }else{
      if(length(file.name)!=nPage){
        warning("Length of 'file.name' is not equal to number of files to be saved, first element in 'file.name' is used",call.=FALSE,immediate.=TRUE)
      }  
      file.name <- file.name[1]   #in case more than one name has been specified

    }
	}else{
      #Each page should be in separate file/window, need multiple file.names
      if(is.null(file.name)){
        file.name <- paste(def.name,"_",ID,sep="")
         #Check that number of file names matches number of pages; if not set file.name to genomeplot 1, .., nPage
        if(length(file.name)!=nPage){
          file.name <- paste(def.name,1:nPage,sep="")
        }
      }else{
        if(length(file.name)<nPage){
          file.name <- paste(file.name[1],1:nPage,sep="")
          warning("Length of 'file.name' is less than number of files to be saved, first element in 'file.name' is reused",call.=FALSE,immediate.=TRUE)
        }else{
          if(length(file.name)>nPage){
            file.name <- file.name[1:nPage]
          }
        }
      }#endif
     
 }#endif
 
 return(file.name)
 
}#endfunction