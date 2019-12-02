####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that finds p-arm and chromosome stopping positions based on cytoband information

##Input:
### cyto.data: dataframe with cytoband information
### unit: the unit used to represent positions in data to be plotted (bp,kbp,mbp)

##Output:
### pstop: a vector giving the stopping positions for each p-arm (adjusted to match unit)
### chromstop: a vector giving the stopping position for each chromosome (adjusted to match unit)

##Required by:
### getArms
### addChromlines
### getGlobPos
### getGlobal.xlim
### adjustSeg

##Requires : none

getArmandChromStop <- function(cyto.data, unit){

	#Sort cyto.data by chromosome number; let be represented by X=23 and Y=24:
	chrom <- cyto.data[,1]
	use.chrom <- gsub("chr","",chrom)  #Remove 'chr' from chrom-strings
	use.chrom[use.chrom=="X"] <- "23"	#Replace X by 23
	use.chrom[use.chrom=="Y"] <- "24"	#Replace Y by 24
	num.chrom <- as.numeric(use.chrom)	#Convert to numeric chromosomes
	
	#Order such that chromosomes are in increasing order from 1:24:
	ord.chrom <- order(num.chrom)
	cyto.data <- cyto.data[ord.chrom,,drop=FALSE] 	

	#Get chromosome stopping positions:
	chrom <- cyto.data[,1]
	chrom.stop <- which(chrom[1:length(chrom)-1]!=chrom[2:length(chrom)])
	chrom.stop <- c(chrom.stop,length(chrom))  #include last chromstop as well

	#Get p-arm stopping positions:
	arm.char <- substring(cyto.data[,4],1,1)   #Retrive first character in name which identifies p and q arms
	arm.stop <- which(arm.char[1:length(arm.char)-1]!=arm.char[2:length(arm.char)])
	p.stop <- arm.stop[-which(arm.stop%in%chrom.stop)]  #Remove qstops

	pos.chromstop <- cyto.data[chrom.stop,3]  #Local stopping position for each chromosome
	pos.pstop <- cyto.data[p.stop,3]		#Local stopping position for each p-arm

	#Factor used to convert positions into desired unit
	f <- switch(unit,
		bp = 1,
		kbp = 10^(-3),
		mbp = 10^(-6))

	return(list(pstop=pos.pstop*f,chromstop=pos.chromstop*f))

}
