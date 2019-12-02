####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################
	
#Function to plot copy numbers and/or pcf segmentation results for one given chromosome with a separate plot for each sample

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
### getPlotParameters
### getFilename
### plotIdeogram
### plotObs
### plotSegments


plotChrom <- function(data=NULL,segments=NULL,pos.unit="bp",sample=NULL,chrom=NULL,assembly="hg19",winsoutliers=NULL,xaxis="pos",layout=c(1,1), plot.ideo=TRUE,...){
	
  #Check, modify and retrieve plot input:
  input <- checkAndRetrievePlotInput(data=data,segments=segments,winsoutliers=winsoutliers,type="chromosome",xaxis=xaxis,pos.unit=pos.unit,sample=sample,chrom=chrom)              
  data <- input$data
  segments <- input$segments
  sampleID <- input$sampleID
  chrom <- input$chrom 
  winsoutliers <- input$winsoutliers
  
	nSeg <- length(segments)   #will be 0 if segments=NULL
	nChrom <- length(chrom)
	nSample <- length(sampleID)
	sample.names <- colnames(data)[-c(1:2)]  #will be NULL if data=NULL
  
	#Plot layout	
	nr <- layout[1]
	nc <- layout[2]	
	
	#If xaxis is index; cannot plot ideogram:
  if(xaxis=="index"){
	 plot.ideo <- FALSE
	}
	
	#Set default plot parameters and change these if user has specified other choices via ... :
	arg <- getPlotParameters(type="chromosome",cr=nc*nr,chrom=chrom,nSeg=nSeg,plot.ideo=plot.ideo,xaxis=xaxis,assembly=assembly,...)
	
	#Check if there should be more than one file/window with plot(s), and get file.name accordingly
	nPage <- ifelse(arg$onefile,1,nChrom)
	file.name <- getFilename(nPage,arg$file.name,ID=paste("Chrom",chrom,sep=" "),type="chromosome")
			
	#Outer argins used for the plot window:
	if(arg$title[1]==""){
		oma <- c(0,0,0,0)
	}else{
		oma <- c(0,0,1,0)
	}			
	mar=c(0.2,0.2,0.3,0.2)
	
	
	#Make it possible to plot more than one chromosome; each chromosomeplot will be in separate file/window:	
	for(c in 1:nChrom){
	  
		#Start new window/file:
		if(!arg$onefile || c==1){
		  #Either print to file, or plot on screen
		  if(!is.null(arg$dir.print)){
			 pdf(file=paste(arg$dir.print,"/",file.name[c],".pdf",sep=""),width=arg$plot.size[1],height=arg$plot.size[2],onefile=TRUE,paper="a4")  #a4-paper 
		  }else{
        if(dev.cur()<=c){       #to make Sweave work
		      dev.new(width=arg$plot.size[1],height=arg$plot.size[2],record=TRUE)
        }
		  }
    }else{
      #Start new page when prompted by user:
		  if(is.null(arg$dir.print)){
		    devAskNewPage(ask = TRUE)   
		  }
    }
		#Initialize
		row=1
		clm=1
		new = FALSE
		
		#Divide the plotting window by the function "framedim":
		frames <- framedim(nr,nc)

    #Which chromosome is this:
		k <- chrom[c]
		ind.chrom <- which(data[,1]==k)  #will be integer(0) if data=NULL
		
		#Get maximum position on chromosome from assembly info
		xmax <- NULL
    if(plot.ideo){
      xmax <- chromMax(chrom=k,cyto.data=arg$assembly,pos.unit=arg$plot.unit)
    }
		
		#Get data ylimits for this chromosome using all sampleID (equalrange=TRUE)
		if(!is.null(data) && arg$equalRange){
      all.sample <- which(sample.names %in% sampleID)
      data.lim <- quantile(data[ind.chrom,all.sample+2],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)	
		}
		#Picking out all segments where chromosome number (in second column) is k, returns new list
		if(!is.null(segments)){
      chrom.segments <- lapply(segments,function(seg,k){seg[seg[,2]==k,]},k=k)    
		}	
		
		#Separate plots for each sample
		
		for(i in 1:nSample){
 			
 			#Since probe positions are the same across the samples, we only want to print warning about probe positons located outside the 
			#chromosome after the last sample has been plotted:
			
			#Select relevant sample
			id <- sampleID[i]
       	
			#Frame dimensions for this plot:
			fig.c <- c(frames$left[clm],frames$right[clm],frames$bot[row],frames$top[row])
			par(fig=fig.c,new=new,oma=oma,mar=mar)
			frame.c <- list(left=frames$left[clm],right=frames$right[clm],bot=frames$bot[row],top=frames$top[row])
	
	   #Divide frame for this sample into a frame for the actual plot and a frame for the ideogram:
      plot.frame <- frame.c
      plot.frame$bot <- frame.c$bot + (frame.c$top-frame.c$bot)*arg$ideo.frac
      ideo.frame <- frame.c
      ideo.frame$top <- plot.frame$bot
      
	    #Get segment limits:
	    seg.lim <- NULL
      if(!is.null(segments)){
		    #Get min and max values in segments to make sure all are shown in plot
        seg.lim <- sapply(chrom.segments,get.seglim,equalRange=arg$equalRange,sampleID=sampleID)  #matrix with limits for each segmentation for this chromosome, min in row 1, max in row2
        seg.lim <- c(min(seg.lim[1,],na.rm=TRUE),max(seg.lim[2,],na.rm=TRUE))   #Get overall min and max over all segments
      }
			
			#Start plotting:
			
			#PLOT IDEOGRAM
      if(plot.ideo){  
  		  par(fig=unlist(ideo.frame),new=new,mar=arg$mar.i)
 			  plotIdeogram(chrom=k,arg$cyto.text,cyto.data=arg$assembly,cex=arg$cex.cytotext,unit=arg$plot.unit)
  			new=TRUE
  		}
  		
			#PLOT DATA POINTS:
			add = FALSE
			if(!is.null(data)){
			 ind.sample <- which(sample.names==id) 
			 if(!arg$equalRange){
				  #Get data limits this chromosome using only this sample (if equalRange=FALSE)
				  data.lim <- quantile(data[ind.chrom,ind.sample+2],probs=c(arg$q/2,(1-arg$q/2)),names=FALSE,type=4,na.rm=TRUE)
	      } 
				
        plotObs(y=data[ind.chrom,ind.sample+2],pos=data[ind.chrom,2],unit=pos.unit,winsoutliers=winsoutliers[ind.chrom,ind.sample+2],type="chromosome",xaxis=xaxis,
            plot.ideo=plot.ideo,k=k,sampleID=id,frame=plot.frame,new=new,op=arg,data.lim=data.lim,seg.lim=seg.lim,xmax=xmax)
				
			  add = TRUE   ##segment plot will be added on top of data plot
			}else{
        data.lim <- NULL
			}
			
			#Plot segments:
			if(!is.null(segments)){
			
        #Plot all segments in list:
        for(s in 1:nSeg){
          use.segments <- chrom.segments[[s]]
          plotSegments(use.segments[use.segments[,1]==id,,drop=FALSE],type="chromosome",k=k,sampleID=id,xaxis=xaxis,add=add,plot.ideo=plot.ideo,col=arg$seg.col[s],
            lty=arg$seg.lty[s],lwd=arg$seg.lwd[s],frame=plot.frame,new=new,unit=pos.unit,seg.lim=seg.lim,data.lim=data.lim,op=arg)
          
          add=TRUE  
        }#endfor
					
        #Add segmentation legends:
        if(!is.null(arg$legend)){
          legend("topright",legend=arg$legend,col=arg$seg.col,lty=arg$seg.lty,cex=arg$cex.axis)
        }  
			
			}
  		
      #If page is full; plot on new page
			if(i%%(nr*nc)==0){
				#Add main title to window page:
				title(arg$title[c],outer=TRUE)
				
				#Start new window page when prompted by user:
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
		title(arg$title[c],outer=TRUE)

		#Close graphcis:
		if(!is.null(arg$dir.print)){
		  if(!arg$onefile){
			 cat("Plot was saved in ",paste(arg$dir.print,"/",file.name[c],".pdf",sep=""),"\n")
			 graphics.off()
      }else{
        if(c==nChrom){
          cat("Plot was saved in ",paste(arg$dir.print,"/",file.name,".pdf",sep=""),"\n") 
          graphics.off()
        }
      }
    }#endif
		
	}#endfor
	
	
}#endfunction






	
	