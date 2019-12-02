
####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to plot observed logR-data, BAF-data and/or segmentation results separately for each sample with chromosomes in different panels

#Input:
### logR: data/frame or matrix with chromosomes in first column, positions in second column and logR-data for one or more samples in subsequent columns
### BAF:  data/frame or matrix with chromosomes in first column, positions in second column and BAF-data for one or more samples in subsequent columns  (must have same dimensions as logR)
### segments: a data frame or a list with data frames containing segmentation results from aspcf
### pos.unit: the unit used to represent the positions in logR
### baf.thres: a vector of length two giving the thresholds above/below which BAF-data are not to be plotted (homozygotes)
### sample: a numeric vector indicating which samples are to be plotted (numerated such that sample=1 indicates the first sample found in logR)
### chrom: numeric/character vector giving the chromosomes to be plotted
### xaxis: what is to be plotted along xaxis; either positions ("pos") or indeces ("index")
### layout: the grid layout for the plot (number of columns and rows)
### plot.ideo: should ideogram be plotted below plots
### ... : other optional plot parameters


##Required by:
### none


##Requires:
### checkAndRetrievePlotInput
### chromMax
### framedim
### get.seglim
### getFilename
### getPlotParameters
### plotIdeogram
### plotObs
### plotSegments

plotAllele <- function(logR=NULL,BAF=NULL,segments=NULL,pos.unit="bp",sample=NULL,chrom=NULL,assembly="hg19",baf.thres=c(0.1,0.9),winsoutliers=NULL,xaxis="pos",layout=c(1,1),plot.ideo=TRUE,...){
  
  #Check plot input, and retrieve desired information:
  data <- list(logR=logR,BAF=BAF)
  input <- checkAndRetrievePlotInput(data=data,segments=segments,winsoutliers=winsoutliers,type="aspcf",xaxis=xaxis,pos.unit=pos.unit,sample=sample,chrom=chrom)
  logR <- input$data$logR
  BAF <- input$data$BAF[,-c(1:2)]   #do not need the first two columns since these are the same as in logR
  segments <- input$segments
  sampleID <- input$sampleID
  chrom <- input$chrom 
  winsoutliers <- input$winsoutliers
 
  
	nSample <- length(sampleID)
	nChrom <- length(chrom)
	nSeg <- length(segments)   #will be 0 if segments=NULL
  sample.names <- colnames(logR)[-c(1:2)]  #will be NULL if data=NULL
  
  #Filter out BAF-values outside thresholds:
  if(!is.null(BAF)){
    BAF[BAF<baf.thres[1]] <- NA	
    BAF[BAF>baf.thres[2]] <- NA
	}
	
  #Plot layout (number of columns and rows in plot grid) specified by user:	
	nr <- layout[1]
	nc <- layout[2]	
		
	#If xaxis is index; cannot plot ideogram:
  if(xaxis=="index"){
	 plot.ideo <- FALSE
	}	
		
	#Set default plot parameters and change these if user has specified other choices via ... :
	arg <- getPlotParameters(type="aspcf",nSeg=nSeg,cr=nc*nr,sampleID=sampleID,plot.ideo=plot.ideo,xaxis=xaxis,assembly=assembly,...)
		
	#Check if there should be more than one file/window with plot(s), and get file.name accordingly
	nPage <- ifelse(arg$onefile,1,nSample)
	file.name <- getFilename(nPage,arg$file.name,ID=sampleID,type="aspcf")
	
	#Divide the plotting window by the function "framedim":
	frames <- framedim(nr,nc)
	
  #Margins used for the plot window:
  if(arg$title[1]!=""){
    oma <- c(0,0,1,0)
  }else{
    oma <- c(0,0,0,0)
  }
  mar=c(0.2,0.2,0.3,0.2)
	
  #Adjust margins separately for logR and BAF-plots:
	mar1 <- arg$mar
	mar1[1] <- 0.3*arg$f
	mar2 <- arg$mar
	mar2[3] <- 0.3*arg$f

  #Two separate arg-lists for logR and BAF
  arg1 <- modifyList(arg,list(xlab="",ylab=arg$ylab[1],h=arg$h[1],mar=mar1))
  arg2 <- modifyList(arg,list(ylab=arg$ylab[2],h=arg$h[2],main="",mar=mar2,connect=FALSE))
  
 
  #Separate plots for each sample:
	
	for(i in 1:nSample){
    
    #Start new window/file:
		if(!arg$onefile || i==1){
		  #Either print to file, or plot on screen
		  if(!is.null(arg$dir.print)){
			 pdf(file=paste(arg$dir.print,"/",file.name[i],".pdf",sep=""),width=arg$plot.size[1],height=arg$plot.size[2],onefile=TRUE,paper="a4")  #a4-paper 
		  }else{
        if(dev.cur()<=i){       #to make Sweave work
		      dev.new(width=arg$plot.size[1],height=arg$plot.size[2],record=TRUE)
        }
		  }
    }else{
      #Start new page when prompted by user:
		  if(is.null(arg$dir.print)){
		    devAskNewPage(ask = TRUE)   
		  }
    }
		
		#Initialize row and column index:
		row=1
		clm=1
		new = FALSE
		
		#Which sample is this:
		id <- sampleID[i]
		ind.sample <- which(sample.names==id)   #will be integer(0) if data=NULL
		
		#Get data limits for this sample if range should include all chromosomes (equalRange=TRUE)
    if(!is.null(logR) && arg$equalRange){
      all.chrom <- which(logR[,1] %in% chrom)
      data.lim1 <- quantile(logR[all.chrom,ind.sample+2],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)
      data.lim2 <- quantile(BAF[all.chrom,ind.sample],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)
		}
		
		#Picking out all segments where sampleid (in first column) is id, returns new list
    if(!is.null(segments)){
      sample.segments <- lapply(segments,function(seg,id){seg[seg[,1]==id,]},id=id)    
		}
		
    #Make separate plots for each chromosome:
		
		for(c in 1:nChrom){ 
		
		  #Frame dimensions for plot c:
			fig.c <- c(frames$left[clm],frames$right[clm],frames$bot[row],frames$top[row])
		  par(fig=fig.c,new=new,oma=oma,mar=mar)
			frame.c <- list(left=frames$left[clm],right=frames$right[clm],bot=frames$bot[row],top=frames$top[row])
			    
      #Insure that plot region is the same for logR and BAF despite different margins:
      fin <- par("fin")        #This figure's width and height
      fig.lines <- fin[2]/0.2  #Number of horizontal lines in this figure
	    logR.frac <- (fig.lines - fig.lines*arg$ideo.frac + (mar1[1]+mar1[3])- (mar2[1]+mar2[3]))/(2*fig.lines)
      baf.frac <- 1-arg$ideo.frac-logR.frac 
       
      #Divide frame for this chromosome into a frame for the logR plot, a frame for the BAF plot and a frame for the ideogram:
      logR.frame <- frame.c 
      logR.frame$bot <-  frame.c$bot + (frame.c$top-frame.c$bot)*(baf.frac+arg$ideo.frac)
      baf.frame <- frame.c
      baf.frame$top <- logR.frame$bot
      baf.frame$bot <- frame.c$bot + (frame.c$top-frame.c$bot)*arg$ideo.frac
      ideo.frame <- frame.c
      ideo.frame$top <- baf.frame$bot
      

			#Select relevant chromosome number
			k <- chrom[c] 			
			
			if(!is.null(segments)){
		    #Get min and max values in logR-segments and BAF-segments to make sure all are shown in plot
        seg.lim1 <- sapply(sample.segments,get.seglim,equalRange=arg$equalRange,k=k)  #matrix with limits for each logR segmentation for this sample, min in row 1, max in row2
        seg.lim1 <- c(min(seg.lim1[1,],na.rm=TRUE),max(seg.lim1[2,],na.rm=TRUE))   #Get overall min and max over all segments
        
        #BAF-segments limits:
        seg.lim2 <- sapply(sample.segments,get.seglim,equalRange=arg$equalRange,k=k,baf=TRUE)  #matrix with limits for each segment, min in row 1, max in row2
        seg.lim2 <- c((1-max(seg.lim2[2,],na.rm=TRUE)),max(seg.lim2[2,],na.rm=TRUE))   #Get overall min and max over all segments; subtract 0.5 in minimum because we mirror around 0.5 later in plotting
      }else{
        seg.lim1 <- NULL
        seg.lim2 <- NULL
      }
      
      #Get maximum position on chromosome from cytoband info
      if(plot.ideo){
        xmax <- chromMax(chrom=k,cyto.data=arg$assembly,pos.unit=arg$plot.unit)
        #PLOT IDEOGRAM
  			figi <- ideo.frame  
  		  par(fig=figi,new=new,mar=arg$mar.i)
 			  plotIdeogram(chrom=k,arg$cyto.text,cyto.data=arg$assembly,cex=arg$cex.cytotext,unit=arg$plot.unit)
  		  new=TRUE
      }else{
        xmax <- NULL
      }
      

      
      #PLOT logR
			add = FALSE
			if(!is.null(logR)){
        ind.chrom <- which(logR[,1]==k)
        if(!arg$equalRange){
          #Get data limits for this sample using just this chromosome (equalRange=FALSE)
				  data.lim1 <- quantile(logR[ind.chrom,ind.sample+2],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)
        } 
				
				#plot logR
        plotObs(y=logR[ind.chrom,ind.sample+2],pos=logR[ind.chrom,2],unit=pos.unit,winsoutliers=winsoutliers[ind.chrom,ind.sample+2],type="aspcf",xaxis=xaxis,
            plot.ideo=TRUE,k=k,frame=logR.frame,new=new,op=arg1,data.lim=data.lim1,seg.lim=seg.lim1,xmax=xmax)
			    	
			  add = TRUE   ##segment plot will be added on top of data plot
			}else{
        data.lim1 <- NULL
			}

      #ADD/PLOT logR SEGMENTS
			if(!is.null(segments)){
        #Plot all logR-segments in list:
        for(s in 1:nSeg){
          use.segments <- sample.segments[[s]]
          use.segments <- use.segments[,c(1:7)]
          plotSegments(use.segments[use.segments[,2]==k,,drop=FALSE],type="aspcf",k=k,xaxis=xaxis,add=add,plot.ideo=TRUE,col=arg$seg.col[s],
            lty=arg$seg.lty[s],lwd=arg$seg.lwd[s],frame=logR.frame,new=new,unit=pos.unit,seg.lim=seg.lim1,data.lim=data.lim1,op=arg1)
          if(k %in% use.segments[,2]){
            add <- TRUE
          }
        }#endfor
					
        #Add segmentation legends:
        if(!is.null(arg$legend)){
          legend("topright",legend=arg$legend,col=arg$seg.col,lty=arg$seg.lty,cex=arg$cex.axis)
        }  
			}
			
			#PLOT BAF:
			add=FALSE
			if(!is.null(BAF)){
        if(!arg$equalRange){
          #Get data limits for this sample using just this chromosome (equalRange=FALSE)
				  data.lim2 <- quantile(BAF[ind.chrom,ind.sample],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)
        } 
				#plot BAF:
        plotObs(BAF[ind.chrom,ind.sample],pos=logR[ind.chrom,2],unit=pos.unit,winsoutliers=NULL,type="aspcf",xaxis=xaxis,k=k,
						plot.ideo=plot.ideo,frame=baf.frame,new=TRUE,op=arg2,xmax=xmax,seg.lim=seg.lim2,data.lim=data.lim2)	
            	
			  add = TRUE   ##segment plot will be added on top of data plot
			}else{
        data.lim2 <- NULL
			}
			
			#ADD/PLOT BAF-SEGMENTS
			if(!is.null(segments)){
        #Plot all BAF-segments in list:
        for(s in 1:nSeg){
          use.segments <- sample.segments[[s]]
          use.segments <- use.segments[,c(1:6,8)]
          plotSegments(use.segments[use.segments[,2]==k,,drop=FALSE],type="aspcf",k=k,xaxis=xaxis,add=add,plot.ideo=TRUE,col=arg$seg.col[s],
            lty=arg$seg.lty[s],lwd=arg$seg.lwd[s],frame=baf.frame,new=TRUE,unit=pos.unit,seg.lim=seg.lim2,data.lim=data.lim2,op=arg2)
              
          #Mirror segments around 0.5, and plot these as well:
          use.segments[,7] <- 1- use.segments[,7]
          plotSegments(use.segments[use.segments[,2]==k,,drop=FALSE],type="aspcf",k=k,xaxis=xaxis,add=TRUE,plot.ideo=TRUE,col=arg$seg.col[s],
            lty=arg$seg.lty[s],lwd=arg$seg.lwd[s],frame=baf.frame,new=TRUE,unit=pos.unit,seg.lim=seg.lim2,data.lim=data.lim2,op=arg2)  
          
          if(k %in% use.segments[,2]){
            add <- TRUE
          }
        }#endfor
					
        #Add segmentation legends:
        if(!is.null(arg$legend)){
          legend("topright",legend=arg$legend,col=arg$seg.col,lty=arg$seg.lty,cex=arg$cex.axis)
        }  
			
			}
			
			
			
			#If page is full; plot on new page
			if(c%%(nr*nc)==0){
				#Add main title to page:
				title(arg$title[i],outer=TRUE)
				
				#Start new page when prompted by user:
				if(is.null(arg$dir.print)){
					devAskNewPage(ask = TRUE)   
				}
				
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

		#Plot sampleid as title
		
		title(arg$title[i],outer=TRUE)
		
		#Close graphcis:
		if(!is.null(arg$dir.print)){
		  if(!arg$onefile){
			 cat("Plot was saved in ",paste(arg$dir.print,"/",file.name[i],".pdf",sep=""),"\n")
			 graphics.off()
      }else{
        if(i==nSample){
          cat("Plot was saved in ",paste(arg$dir.print,"/",file.name,".pdf",sep=""),"\n") 
          graphics.off()
        }
      }
    }#endif
		
	}#endfor
		
	
}#endfunction



