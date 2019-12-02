####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function that plots the frequency of deletions and amplifications given a threshold - by genome og chromosomes

##Input:
### segments: result from segmentation (could be a data frame with the estimates, or the segments data frame). Also possible to specify a data frame with original data
### thres.gain,thres.loss: a vector with threshold(s) to be applied for abberration calling
### pos.unit: unit used for positions (bp,kbp,mbp)
### chrom: a vector with chromosomes to be plotted. If specified the frequencies are plotted with one panel for each chromosome
### layout: number of columns and rows in plot
### ... : other optional plot parameters


##Required by: none

##Requires: 
### numericChrom
### checkChrom
### genomeFreq
### chromosomeFreq
### pullOutContent
### getFreqData

plotFreq <- function(segments,thres.gain,thres.loss=-thres.gain,pos.unit="bp",chrom=NULL,layout=c(1,1),...){
	
	#Check data input:
	#Make sure data is a data frame (could be a list if it contains segmentation results or results from winsorize)
  data = segments
  if(is.data.frame(data)){
    #Check if segments or data:
    data = getFreqData(data)  
  }else{
    data <- pullOutContent(res=data,what="estimates")
	}
	
  stopifnot(ncol(data)>=3)  #something is missing in input data
	
	#Make sure that chromosomes in data are numeric:
  data[,1] <- numericChrom(data[,1])
  
  #If chrom is unspecified, the whole genome is plotted. Otherwise, selected chromosomes are plotted
  type <- ifelse(is.null(chrom),"genome","bychrom")
  
	#Check input chrom, alt. get unique chromosomes 
	chrom <- checkChrom(data=data,segments=NULL,chrom=chrom)
	
	#Check pos.unit input:
  if(!pos.unit %in% c("bp","kbp","mbp")){
    stop("pos.unit must be one of bp, kbp and mbp",call.=FALSE)
  }
  
  #Making sure number of gain thresholds and loss thresholds are the same:
  nT <- min(length(thres.gain),length(thres.loss))
  thres.gain <- thres.gain[1:nT]
  thres.loss <- thres.loss[1:nT]
 
	#plot frequency, either over genome or by chromosomes:
	switch(type,
		genome = genomeFreq(data,thres.gain,thres.loss,pos.unit,layout,...),
		bychrom = chromosomeFreq(data,thres.gain,thres.loss,pos.unit,chrom,layout,...)
		)

}


#Function that plots frequencies by genome 

##Required by: plotFreq

##Requires: 
### getFreqPlotParameters
### getGlobal.xlim 
### adjustPos
### framedim
### updateFreqParameters
### addToFreqPlot
### addChromlines
### addArmlines
### chromPattern

genomeFreq <- function(data,thres.gain,thres.loss,pos.unit,layout,...){

	nr <- layout[1]
	nc <- layout[2]
	cr <- nr*nc

  nProbe <- nrow(data)
  nT <- length(thres.gain)
  
  
	#Get plot parameters:
	op <- getFreqPlotParameters(type="genome",nc=nc,nr=nr,thres.gain=thres.gain,thres.loss=thres.loss,...)

  #Set global xlimits if not specified by user:
  if(is.null(op$xlim)){
    op$xlim <- getGlobal.xlim(op=op,pos.unit=pos.unit,chrom=unique(data[,1]))
  }
  
  #Adjust positions to be plotted along xaxis; i.e. get global positions, scale according to plot.unit, and get left and right pos for freq.rectangles
  #to be plotted (either continuous or 1 probe long):
  x <- adjustPos(position=data[,2],chromosomes=data[,1],pos.unit=pos.unit,type="genome",op=op)
  xleft <- x$xleft
  xright <- x$xright
  
  
  if(dev.cur()<=1){       #to make Sweave work
    dev.new(width=op$plot.size[1],height=op$plot.size[2],record=TRUE)
  }
	  

  #Initialize:
  row=1
  clm=1
  new = FALSE

  #Division of plotting window:
  frames <- framedim(nr,nc)

  #One plot for each value in thres.gain/thres.loss:
  for(t in 1:nT){

    #Frame dimensions for plot t:
		fig.t <- c(frames$left[clm],frames$right[clm],frames$bot[row],frames$top[row])
		par(fig=fig.t,new=new,oma=c(0,0,1,0),mar=op$mar)
		
    #Calculate the percentage of samples that have estimated copy number larger than thres at given position:
		freq.amp <- rowMeans(data[,-c(1:2),drop=FALSE] > thres.gain[t])*100
		freq.del <- rowMeans(data[,-c(1:2),drop=FALSE] < thres.loss[t])*100

    #Find default ylimits and at.y (tickmarks):
    op <- updateFreqParameters(freq.del,freq.amp,op)
    
		#Empty plot with correct limits
		plot(1,1,type="n",ylim=op$ylim,xlim=op$xlim,xaxs="i",main="",frame.plot=TRUE,yaxt="n",xaxt="n",ylab="",xlab="")

    #Add shifting white/grey pattern to backgroud to separate chromosomes:
    chromPattern(pos.unit,op)
    
    #main title for this plot
    title(main=op$main[t],line=op$main.line,cex.main=op$cex.main)

		#Add axes, labels and percentage lines:
		addToFreqPlot(op,type="genome")

		#Plot frequencies:
    rect(xleft=xleft,ybottom=0,xright=xright,ytop=freq.amp,col=op$col.gain,border=op$col.gain)
    rect(xleft=xleft,ybottom=0,xright=xright,ytop=-freq.del,col=op$col.loss,border=op$col.loss)

    #Add line at y=0:
		abline(h=0,lty=1,col="grey82",lwd=1.5)

    #Add line at y=0:
		abline(h=0,lty=1,col="grey82",lwd=1.5)

		#Separate chromosomes by vertical lines:
		op$chrom.lty = 1
		addChromlines(data[,1],xaxis="pos",unit=pos.unit,cex=op$cex.chrom,op=op)
		addArmlines(data[,1],xaxis="pos",unit=pos.unit,cex=op$cex.chrom,op=op)

		#Get new page, or update column/row:
		if(t%%(nr*nc)==0){
			#Start new plot page (prompted by user)
			devAskNewPage(ask = TRUE)

			#Reset columns and row in layout:
			clm = 1
			row = 1
			new=FALSE
		}else{
			#Update column and row index:
			if(clm<nc){
				clm <- clm+1
			}else{
				clm <- 1
				row <- row+1
			}#endif
			new=TRUE
		}#endif

	}#endfor


}#end function




#Function that plots frequencies by chromosomes (each chrom in separate panel)


##Required by: plotFreq

##Requires: 
### getFreqPlotParameters 
### adjustPos
### framedim
### updateFreqParameters
### plotIdeogram
### chromMax
### addToFreqPlot


chromosomeFreq <- function(data,thres.gain,thres.loss,pos.unit,chrom,layout,...){

  nProbe <- nrow(data)
	nT <- length(thres.gain)	
  nChrom <- length(chrom)

  #Grid layout
	nr <- layout[1]
	nc <- layout[2]

	#Get plot parameters:
	op <- getFreqPlotParameters(type="bychrom",nc=nc,nr=nr,thres.gain=thres.gain,thres.loss=thres.loss,chrom=chrom,...)

	#Margins for entire plot in window:
	if(all(op$title=="")){
		oma <- c(0,0,0,0)
	}else{
		oma <- c(0,0,1,0)
	}
  mar=c(0.2,0.2,0.3,0.2)
  

  #Adjust positions to be plotted along xaxis; i.e. scale according to plot.unit, and get left and right pos for freq.rectangles to be plotted (either continuous
  #or 1 probe long):
  x <- adjustPos(position=data[,2],chromosomes=data[,1],pos.unit=pos.unit,type="chromosome",op=op)
  xleft <- x$xleft
  xright <- x$xright
  
	#Divide the plotting window by the function "framedim":
	frames <- framedim(nr,nc)

	#make separate plots for each value of thres.gain/thres.loss
	for(t in 1:nT){
        
    if(dev.cur()<=1){       #to make Sweave work
      dev.new(width=op$plot.size[1],height=op$plot.size[2],record=TRUE)
    }
		  
    #Initialize row and column index:
    row=1
    clm=1
    new = FALSE

    #For each probe, get percentage of samples amplified and deleted for this thres:
		freq.amp <- rowMeans(data[,-c(1:2),drop=FALSE] > thres.gain[t])*100
		freq.del <- rowMeans(data[,-c(1:2),drop=FALSE] < thres.loss[t])*100
		
		#Find default ylimits and at.y (tickmarks):
		use <- which(data[,1] %in% chrom)
    op <- updateFreqParameters(freq.del[use],freq.amp[use],op)
      
    #Make separate plots for each chromosome:

		for(c in 1:nChrom){

			#Frame dimensions for plot c:
			fig.c <- c(frames$left[clm],frames$right[clm],frames$bot[row],frames$top[row])
			par(fig=fig.c,new=new,oma=oma,mar=mar)
			
			#Make list with frame dimensions:
			frame.c <- list(left=frames$left[clm],right=frames$right[clm],bot=frames$bot[row],top=frames$top[row])

			#Select relevant chromosome number
			k <- chrom[c]

			#Pick out frequencies for this chromsome
			ind.c <- which(data[,1]==k)
			freqamp.c <- freq.amp[ind.c]
			freqdel.c <- freq.del[ind.c]

      xlim <- c(0,max(xright[ind.c])) 
      
      #Plot ideogram at below frequencies:
			if(op$plot.ideo){ 
				#Ideogram frame:
				ideo.frame <- frame.c
				ideo.frame$top <- frame.c$bot + (frame.c$top-frame.c$bot)*op$ideo.frac
		    
        par(fig=unlist(ideo.frame),new=new,mar=op$mar.i)
        #Plot ideogram and get maximum probe position in ideogram:
				plotIdeogram(chrom=k,cyto.text=op$cyto.text,cyto.data=op$assembly,cex=op$cex.cytotext,unit=op$plot.unit)

        #Get maximum position for this chromosome:
        xmaxI <- chromMax(chrom=k,cyto.data=op$assembly,pos.unit=op$plot.unit)
				xlim <- c(0,xmaxI)

				new <- TRUE
			}

			#Freq.plot-dimensions:
			frame.c$bot <- frame.c$bot + (frame.c$top-frame.c$bot)*op$ideo.frac
			par(fig=unlist(frame.c),new=new,mar=op$mar)

      #Limits:
      if(!is.null(op$xlim)){
        xlim <- op$xlim
      }
        
			#Empty plot:
			plot(1,1,type="n",ylim=op$ylim,xlim=xlim,xaxs="i",main=op$main[c],frame.plot=FALSE,yaxt="n",xaxt="n",cex.main=op$cex.main,ylab="",xlab="")

			#Add axes, labels and percentageLines:
			chrom.op <- op
		  chrom.op$xlim <- xlim
			
			addToFreqPlot(chrom.op,type="bychrom")

			#Plot frequencies as rectangles
      rect(xleft=xleft[ind.c],ybottom=0,xright=xright[ind.c],ytop=freqamp.c,col=op$col.gain,border=op$col.gain) 
      rect(xleft=xleft[ind.c],ybottom=0,xright=xright[ind.c],ytop=-freqdel.c,col=op$col.loss,border=op$col.loss)
      

			#Add line at y=0 and x=0
			abline(h=0,lty=1,col="grey90")
			abline(v=0)



			#If page is full; start plotting on new page
			if(c%%(nr*nc)==0 && c!=nChrom){
				#Add main title to page:
				title(op$title[t],outer=TRUE)

				#Start new page when prompted by user:
				devAskNewPage(ask = TRUE)   

				#Reset columns and row in layout:
				clm = 1
				row = 1
				new=FALSE

			}else{
				#Update column and row index:
				if(clm<nc){
					clm <- clm+1
				}else{
					clm <- 1
					row <- row+1
				}#endif
				new=TRUE
			}#endif


		}#endfor

    #Add main title to page:
    title(op$title[t],outer=TRUE)
   	if(t!=nT){
		  devAskNewPage(ask = TRUE)
		}
		
	}#endfor


}#endfunction

		







