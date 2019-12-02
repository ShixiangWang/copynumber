####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################

#Function to plot observed copy number data and/or segmentation result for one given sample with a separate figure for each chromosome

#Input:
### data: dataframe or matrix with chromosomes in first column, positions in second column and copy number data for one or more samples in subsequent columns
### segments: a data frame or a list with data frames containing segmentation results
### pos.unit: the unit used to represent the positions in data
### sample: a numeric vector indicating which samples are to be plotted (numerated such that sample=1 indicates the first sample found in data)
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


plotSample <- function(data=NULL,segments=NULL,pos.unit="bp",sample=NULL,chrom=NULL,assembly="hg19",winsoutliers=NULL,xaxis="pos",layout=c(1,1),plot.ideo=TRUE,...){
                
  #Check, modify and retrieve plot input:
  input <- checkAndRetrievePlotInput(data=data,segments=segments,winsoutliers=winsoutliers,type="sample",xaxis=xaxis,pos.unit=pos.unit,sample=sample,chrom=chrom)              
  data <- input$data
  segments <- input$segments
  sampleID <- input$sampleID
  chrom <- input$chrom
  winsoutliers <- input$winsoutliers
  
  
	nSample <- length(sampleID)
	nChrom <- length(chrom)
	nSeg <- length(segments)   #will be 0 if segments=NULL
  sample.names <- colnames(data)[-c(1:2)]  #will be NULL if data=NULL
  	
	#Plot layout (number of columns and rows in plot grid) specified by user:	
	nr <- layout[1]
	nc <- layout[2]	
		
	#If xaxis is index; cannot plot ideogram:
  if(xaxis=="index"){
	 plot.ideo <- FALSE
	}	
	#Set default plot parameters and change these if user has specified other choices via ... :
	arg <- getPlotParameters(type="sample",nSeg=nSeg,cr=nc*nr,sampleID=sampleID,plot.ideo=plot.ideo,xaxis=xaxis,assembly=assembly,...)
		
	#Margins used for the plot window:
	if(arg$title[1]==""){
		oma <- c(0,0,0,0)
	}else{
		oma <- c(0,0,1,0)
	}
	mar=c(0.2,0.2,0.3,0.2)
		
	
	#Check if there should be more than one file/window with plot(s), and get file.name accordingly
	nPage <- ifelse(arg$onefile,1,nSample)
	file.name <- getFilename(nPage,arg$file.name,ID=sampleID,type="sample")
	
	#Divide the plotting window by the function "framedim":
	frames <- framedim(nr,nc)
		
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
    if(!is.null(data) && arg$equalRange){
      all.chrom <- which(data[,1] %in% chrom)
      data.lim <- quantile(data[all.chrom,ind.sample+2],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)
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
	
	    #Divide frame for this chromosome into a frame for the actual plot and a frame for the ideogram:
      plot.frame <- frame.c
      plot.frame$bot <- frame.c$bot + (frame.c$top-frame.c$bot)*arg$ideo.frac
      ideo.frame <- frame.c
      ideo.frame$top <- plot.frame$bot
      
			#Select relevant chromosome number
			k <- chrom[c] 			
			
			if(!is.null(segments)){
		    #Get min and max values in segments to make sure all are shown in plot
        seg.lim <- sapply(sample.segments,get.seglim,equalRange=arg$equalRange,k=k)  #matrix with limits for each segmentation for this sample, min in row 1, max in row2
        seg.lim <- c(min(seg.lim[1,],na.rm=TRUE),max(seg.lim[2,],na.rm=TRUE))   #Get overall min and max over all segments
      }else{
        seg.lim <- NULL
      }
      
      #Get maximum position on chromosome from cytoband info
      if(plot.ideo){
        xmax <- chromMax(chrom=k,cyto.data=arg$assembly,pos.unit=arg$plot.unit)
        #PLOT IDEOGRAM  
		    par(fig=unlist(ideo.frame),new=new,mar=arg$mar.i)
		    plotIdeogram(chrom=k,arg$cyto.text,cyto.data=arg$assembly,cex=arg$cex.cytotext,unit=arg$plot.unit)
		    new=TRUE
      }else{
        xmax <- NULL
      }
      
			#PLOT DATA POINTS:
			add = FALSE                                                                                           
			if(!is.null(data)){
			 ind.chrom <- which(data[,1]==k)
			 if(!arg$equalRange){
				  #Get data limits for this sample using just this chromosome (equalRange=FALSE)
				  data.lim <- quantile(data[ind.chrom,ind.sample+2],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)
	      } 
				
        plotObs(y=data[ind.chrom,ind.sample+2],pos=data[ind.chrom,2],unit=pos.unit,winsoutliers=winsoutliers[ind.chrom,ind.sample+2],type="sample",xaxis=xaxis,
            plot.ideo=plot.ideo,k=k,frame=plot.frame,new=new,op=arg,data.lim=data.lim,seg.lim=seg.lim,xmax=xmax)
				
			  add = TRUE   ##segment plot will be added on top of data plot
			}else{
        data.lim <- NULL
			}
			
			#Plot segments:
			if(!is.null(segments)){
			
        #Plot all segments in list:
        for(s in 1:nSeg){
          use.segments <- sample.segments[[s]]
          plotSegments(use.segments[use.segments[,2]==k,,drop=FALSE],type="sample",k=k,xaxis=xaxis,add=add,plot.ideo=plot.ideo,col=arg$seg.col[s],
            lty=arg$seg.lty[s],lwd=arg$seg.lwd[s],frame=plot.frame,new=new,unit=pos.unit,seg.lim=seg.lim,data.lim=data.lim,op=arg)
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



