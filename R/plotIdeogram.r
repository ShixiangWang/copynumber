####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that plots ideogram for given chromosome

#Input:
### chrom: number between 1-24 indicating which chromosome's ideogram is to be plotted
### cyto.text: should cytoband-names be plotted along ideogram?
### cex: the size used for cyto.text
### cyto.data: data frame with cytoband-information
### cyto.unit: the unit used to represent positons in cyto.data
### unit: the unit used for positions in the plot


##Required by:
### plotFreq (chromosomeFreq)
### plotAllele
### plotChrom
### plotSample
### plotHeatmap
### plotWeightedFreq

##Requires:
### convert.unit

plotIdeogram <- function(chrom,cyto.text=FALSE,cex=0.6,cyto.data,cyto.unit="bp",unit){
	
	if(chrom==23){
		chrom.cytoband <- cyto.data[cyto.data[,1]=="chrX",]
	}else{
		if(chrom==24){
			chrom.cytoband <- cyto.data[cyto.data[,1]=="chrY",]
		}else{
			chrom.cytoband <- cyto.data[cyto.data[,1]==paste("chr",chrom,sep=""),]
		}
	}
	
	cyto.start <- chrom.cytoband[,2]
	cyto.end <- chrom.cytoband[,3]
	scale <- convert.unit(unit1=unit,unit2=cyto.unit)
	
	xleft <- cyto.start*scale
	xright <- cyto.end*scale
	n <- length(xleft)
	chrom.length <- xright[n]-xleft[1]
	
	stain <- chrom.cytoband[,5]
	sep.stain <- c("gpos","gneg","acen","gvar","stalk")

	g <- sapply(sep.stain,grep,x=stain,fixed=TRUE)

	centromere <- g$acen
	stalk <- g$stalk
	col <- rep("",n)
	col[stain=="gneg"] <- "white"
	col[stain=="gpos100"] <- "black"
	col[stain=="gpos75"] <- "gray25"
	col[stain=="gpos50"] <- "gray50"
	col[stain=="gpos25"] <- "gray75"
	col[stain=="stalk"] <- "gray90"
	col[stain=="gvar"] <- "grey"
	col[stain=="acen"] <- "yellow"
	density <- rep(NA,n)
	angle <- rep(45,n)
	density[stain=="gvar"] <- 15
	

	ylow <- 0
	yhigh <- 1
	

	plot(x=c(0,max(xright)),y=c(ylow,yhigh),type="n",axes=FALSE,xlab="",ylab="",xlim=c(0,max(xright)),ylim=c(0,1),xaxs="i")
	
	#Rectangles:
	skip.rect <- c(1,centromere,n,stalk)
	rect(xleft[-skip.rect],rep(ylow,n-length(skip.rect)),xright[-skip.rect],rep(yhigh,n-length(skip.rect)),
		col=col[-skip.rect],border="black",density=density[-skip.rect],angle=angle[-skip.rect])
	
	#Round edges at ideogram start, stop and at centromere:
	draw.roundEdge(start=xleft[1],stop=xright[1],y0=ylow,y1=yhigh,col=col[1],bow="left",density=density[1],angle=angle[1],chrom.length=chrom.length)
	draw.roundEdge(start=xleft[centromere[1]],stop=xright[centromere[1]],y0=ylow,y1=yhigh,col=col[centromere[1]],bow="right",density=density[centromere[1]],
		angle=angle[centromere[1]],lwd=1,chrom.length=chrom.length)
	draw.roundEdge(start=xleft[centromere[2]],stop=xright[centromere[2]],y0=ylow,y1=yhigh,col=col[centromere[2]],bow="left",density=density[centromere[2]],
		angle=angle[centromere[2]],lwd=1,chrom.length=chrom.length)
	draw.roundEdge(start=xleft[n],stop=xright[n],y0=ylow,y1=yhigh,col=col[n],bow="right",density=density[n],angle=angle[n],chrom.length=chrom.length)

	#Draw stalk-segment:
	if(length(stalk)>0){
		for(i in 1:length(stalk)){
			drawStalk(xleft[stalk[i]],xright[stalk[i]],ylow,yhigh,col=col[stalk[i]])
		}
	}
	if(cyto.text){
		mtext(text=paste(chrom.cytoband[,4],"-",sep=" "),side=1,at=(xleft + (xright-xleft)/2),cex=cex,las=2,adj=1,xpd=NA)#,line=-1)#,outer=TRUE)
	}
	
}




draw.roundEdge <- function(start,stop,y0,y1,col,bow,density=NA,angle=45,lwd=1,chrom.length){
	#Y points in round edge:
	f <- rep(0,0)
	f[1] <- 0.001
	i=1
	half <- y0+(y1-y0)/2
	while(f[i]<half){
		f[i+1] <- f[i]*1.3
		i <- i+1
	}
	f <- f[-length(f)]
	
	Y <- c(y1,y1,y1-f,half,y0+rev(f),y0,y0)
	
	#X points in roundedge
	cyto.length <- stop-start
	
	share <- cyto.length/chrom.length
	if(share>0.2){
		#to create bow in end of chromosome 24
		share <- 0.2
	}
	
	if(bow=="left"){
	
		round.start <- start + cyto.length*(1-share)^20
		
		x <- seq(round.start,start,length.out=(length(f)+2))
		revx <- rev(x[-length(x)])
		x <- c(x,revx)
		X <- c(stop,x,stop)
	}else{
		if(bow=="right"){
			round.start <- stop - cyto.length*(1-share)^20
			x <- seq(round.start,stop,length.out=(length(f)+2))
			revx <- rev(x[-length(x)])
			x <- c(x,revx)
			X <- c(start,x,start)

		}
	}
	
	polygon(x=X,y=Y,col=col,border="black",density=density,angle=angle,lwd=lwd)

}

drawStalk <- function(start,stop,y0,y1,col){
	size <- stop-start
	x1 <- c(start,start+size/3,stop-size/3,stop)
	x2 <- rev(x1)
	x <- c(x1,x2)
	y_1 <- c(y0,y0+0.25,y0+0.25,y0)
	y_2 <- c(y1,y1-0.25,y1-0.25,y1)
	y <- c(y_1,y_2)
	polygon(x=x,y=y,col=col)

}