####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to plot observed copy number estimates and/or segmentation results across entire genome

#Input:
### data: dataframe or matrix with chromosomes in first column, positions in second column and copy number data for one or more samples in subsequent columns
### segments: a data frame or a list with data frames containing segmentation results
### pos.unit: the unit used to represent the positions in data
### sample: a numeric vector indicating which samples are to be plotted (numerated such that sample=1 indicates the first sample found in data)
### winsoutliers: a dataframe/matrix of same dimensions as data indicating outliers statuses
### xaxis: what is to be plotted along xaxis; either positions ("pos") or indeces ("index")
### layout: the grid layout for the plot (number of columns and rows)
### ... : other optional plot parameters


##Required by:
### none

##Requires:
### checkAndRetrievePlotInput
### framedim
### get.seglim
### getFilename
### getGlobal.xlim
### getPlotParameters
### plotObs
### plotSegments


plotGenome <- function(data=NULL,segments=NULL,pos.unit="bp",sample=NULL,assembly="hg19",winsoutliers=NULL,xaxis="pos",layout=c(1,1),...){
	
	#Check, modify and retrieve plot input:
  input <- checkAndRetrievePlotInput(data=data,segments=segments,winsoutliers=winsoutliers,type="genome",xaxis=xaxis,pos.unit=pos.unit,sample=sample)              
  data <- input$data
  segments <- input$segments
  sampleID <- input$sampleID
  winsoutliers <- input$winsoutliers
 
		
	nSample <- length(sampleID)	
	nSeg <- length(segments)   #will be 0 if segments=NULL
  sample.names <- colnames(data)[-c(1:2)]  #will be NULL if data=NULL
  	
    	
	#Plot layout	
	nr <- layout[1]
	nc <- layout[2]	
	rc <- nr*nc		
      	
	#Set default plot parameters and change these if user has specified other choices via ... :
	arg <- getPlotParameters(type="genome",cr=nc*nr,nSeg=nSeg,sampleID=sampleID,plot.ideo=FALSE,xaxis=xaxis,assembly=assembly,...)   
	
 #Set global xlimits if not specified by user:
  if(is.null(arg$xlim) && xaxis=="pos"){
    if(!is.null(data)){
      chr <- unique(data[,1])
    }else{
      chr <- unique(unlist(sapply(segments,function(seg){unique(seg[,2])})))
    }
    arg$xlim <- getGlobal.xlim(op=arg,pos.unit=pos.unit,chrom=chr)
  }
  
 #Get data limits if equalRange -> range will be max and min across all samples
  if(!is.null(data) && arg$equalRange){
    all.sample <- which(sample.names %in% sampleID)
    data.lim <- quantile(data[,all.sample+2],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)
	}	
 
	#Check if there should be more than one file/window with plot(s), and get file.name accordingly
	nPage <- ifelse(arg$onefile,1,ceiling(nSample/(rc)))
	file.name <- getFilename(nPage,arg$file.name,ID=sampleID,type="genome")
	
  #New window/pdf-file for each file.name
  nPlot <- 1
  	
  for(j in 1:length(file.name)){
    #Either print to file, or plot on screen		
    if(!is.null(arg$dir.print)){
		  pdf(file=paste(arg$dir.print,"/",file.name[j],".pdf",sep=""),width=arg$plot.size[1],height=arg$plot.size[2],onefile=TRUE,paper="a4")  #a4-paper 
	  }else{
		  #windows(width=arg$plot.size[1],height=arg$plot.size[2],record=TRUE)
		  if(dev.cur()<=j){       #to make Sweave work
		    dev.new(width=arg$plot.size[1],height=arg$plot.size[2],record=TRUE)
      }
   }

    #Initialize:
	  row=1
    clm=1
    new = FALSE

    #Division of plotting window:
    frames <- framedim(nr,nc)

    if(!arg$onefile){
      #Pick out samples that are to be plotted in this window/file:
      use.sampleID <- sampleID[nPlot:(rc*j)]  
      use.sampleID <- use.sampleID[!is.na(use.sampleID)]  #if number of samples divided by rc is not an integer we will get NA's here at some point
    }else{
      use.sampleID <- sampleID
    }	
    
    
    #Separate plots for each sample
  
  
    for(i in 1:length(use.sampleID)){
			
			#Frame dimensions for plot i:
			fig.c <- c(frames$left[clm],frames$right[clm],frames$bot[row],frames$top[row])
			par(fig=fig.c,new=new,oma=c(0,0,0.5,0),mar=arg$mar)
			frame.c <- list(left=frames$left[clm],right=frames$right[clm],bot=frames$bot[row],top=frames$top[row])

      #Which sample is this:
		  id <- use.sampleID[i]
		  ind.sample <- which(sample.names==id)   #will be integer(0) if data=NULL
		  
		  #Get min and max values in segments to make sure all are shown in plot
      if(!is.null(segments)){
        seg.lim <- sapply(segments,get.seglim,equalRange=arg$equalRange,sampleID=id)  #matrix with limits for each segment, min in row 1, max in row2
        seg.lim <- c(min(seg.lim[1,]),max(seg.lim[2,]))   #Get overall min and max over all segments         
	    }else{
        seg.lim <- NULL
      }
	    
      
      #PLOT DATA POINTS:
			add = FALSE
			if(!is.null(data)){
        if(!arg$equalRange){
				  #Get data limits for just this sample
				  data.lim <- quantile(data[,ind.sample+2],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)
        } 
				                                 
        plotObs(y=data[,ind.sample+2],pos=data[,2],unit=pos.unit,winsoutliers=winsoutliers[,ind.sample+2],type="genome",xaxis=xaxis,sampleID=id,
            chromosomes=data[,1],frame=frame.c,new=new,op=arg,data.lim=data.lim,seg.lim=seg.lim)
				
			  add = TRUE   ##segment plot will be added on top of data plot
			}else{
        data.lim <- NULL
			}
			
			#PLOT SEGMENTS:
			if(!is.null(segments)){
			  sample.segments <- lapply(segments,function(seg,id){seg[seg[,1]==id,]},id=id)    
        #Plot all segments in list:
        for(s in 1:nSeg){
          use.segments <- sample.segments[[s]]
          plotSegments(use.segments,type="genome",xaxis=xaxis,add=add,col=arg$seg.col[s],sampleID=id,
            lty=arg$seg.lty[s],lwd=arg$seg.lwd[s],frame=frame.c,new=new,unit=pos.unit,seg.lim=seg.lim,data.lim=data.lim,op=arg)
          
          add <- TRUE
        }#endfor
					
        #Add segmentation legends:
        if(!is.null(arg$legend)){
          legend("topright",legend=arg$legend,col=arg$seg.col,lty=arg$seg.lty,cex=arg$cex.axis)
        }  
			
			}
		
			if(i%%(nr*nc)==0){
			  #Add main title to window page:
				title(arg$title,outer=TRUE)
				
				#Start new plot page (prompted by user)
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
	
	  #Add main title to page:
		title(arg$title,outer=TRUE)
		
		#Close graphics
		if(!is.null(arg$dir.print)){
			cat("Plot was saved in ",paste(arg$dir.print,"/",file.name[j],".pdf",sep=""),"\n")
			graphics.off()
		}
		
		#Update current number of plots:
		nPlot <- rc*j +1 
		
  }#endfor

	
}#endfunction







