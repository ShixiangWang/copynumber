
####################################################################
## Author: Gro Nilsen, Knut Liestøl and Ole Christian Lingjærde.
## Maintainer: Gro Nilsen <gronilse@ifi.uio.no>
## License: Artistic 2.0
## Part of the copynumber package
## Reference: Nilsen and Liestøl et al. (2012), BMC Genomics
####################################################################


##Required by:


##Requires:
### findNN
### getMad
### pcf(not the main function, but the rest of the help functions found in the same document)
### handleMissing
### pullOutContent


## Main function for pcf-analysis to be called by the user

pcfPlain <- function(pos.data,kmin=5,gamma=40,normalize=TRUE,fast=TRUE,digits=4,return.est=FALSE,verbose=TRUE){
  
  
  #Input could come from winsorize and thus be a list; check and possibly retrieve data frame wins.data
  pos.data <- pullOutContent(pos.data,what="wins.data")
    
  #Make sure all data columns are numeric:
  if(any(!sapply(pos.data,is.numeric))){
    stop("All input in pos.data must be numeric",call.=FALSE)
  }
  
  #Check data input:
  stopifnot(ncol(pos.data)>=2)  #something is missing in input data
  
  #Extract information from pos.data:
  position <- pos.data[,1]
  cn.data <- as.matrix(pos.data[,-1])
  nSample <- ncol(pos.data)-1
  sampleid <- colnames(pos.data)[-1]  
  nProbe <- length(position)
	
  
  #save user's gamma
	gamma0 <- gamma     
	
	sd <- rep(1,nSample) #sd is only used if normalize=TRUE, and then these values are replaced by MAD-sd
	#If number of probes in entire data set is less than 100K, the MAD sd-estimate is calculated using all obs for every sample
  #Only required if normalize=T
  if(nProbe<100000 && normalize){
    #calculate MAD-sd for each sample:
    for(j in 1:nSample){
      sample.data <- pos.data[,j+1]
      sd[j] <- getMad(sample.data[!is.na(sample.data)],k=25)   #Take out missing values before calculating mad
		}
  }#endif
  
  
	#Initialize
	pcf.names <- c("pos",sampleid)
	seg.names <- c("sampleID","start.pos","end.pos","n.probes","mean")
	
  segments <- data.frame(matrix(nrow=0,ncol=5))
  colnames(segments) <- seg.names
  if(return.est){
    pcf.est <- matrix(nrow=0,ncol=nSample)
	}
  
  		
  #Run PCF separately for each sample:
  for(i in 1:nSample){
    if(return.est){
      #Initialize:
      yhat <- rep(NA,length(nProbe))
    }
    
    sample.data <- cn.data[,i]
    
		#Remove probes with missing obs; Only run pcf on non-missing values
		obs <- !is.na(sample.data)
		obs.data <- sample.data[obs]
	
	  
	  if(length(obs.data)>0){  ##Make sure there are observations for this sample! If not, estimates are left NA as well
		
			#If number of probes in entire data set is >= 100K, the MAD sd-estimate is calculated using obs in this arm for this sample.
      #Only required if normalize=T 
      if(nProbe>=100000 && normalize){
        sd[i] <- getMad(obs.data,k=25)   
      }
    
      #Scale gamma by variance if normalize is TRUE 
      use.gamma <- gamma0
      if(normalize){
        use.gamma <- gamma0*(sd[i])^2
      }
      
  		#Must check that sd!=0 and sd!!=NA -> use.gamma=0/NA. If not, simply calculate mean of observations
			if(use.gamma==0 || is.na(use.gamma)){  
			  if(return.est){
          res <- list(Lengde=length(obs.data),sta=1,mean=mean(obs.data),nIntervals=1,yhat=rep(mean(obs.data)))
        }else{
          res <- list(Lengde=length(obs.data),sta=1,mean=mean(obs.data),nIntervals=1)
        }
        
			}else{
  			# Compute piecewise constant fit
  			#run fast approximate PCF if fast=TRUE and number of probes>400, or exact PCF otherwise
  			if(!fast || length(obs.data)<400){
  				#Exact PCF:	
  			 	res <- exactPcf(y=obs.data,kmin=kmin,gamma=use.gamma,yest=return.est)	 			
  			}else{
          #Run fast PCF:
          res <- selectFastPcf(x=obs.data,kmin=kmin,gamma=use.gamma,yest=return.est)
  			}#endif
  			
  			
			}#endif
			
			#Retrieve segment info from results:
			seg.start <- res$sta
		  seg.stop <- c(seg.start[-1]-1,length(obs.data))
		  seg.npos <- res$Lengde
			seg.mean <- res$mean
			nSeg <- res$nIntervals 
      if(return.est){
        yhat[obs] <- res$yhat
      } 
      #Find genomic positions for start and stop of each segment:
		  pos.start <- position[seg.start]
		  pos.stop <- position[seg.stop]
			
      #Handle missing values:
			if(any(!obs)){
			  #first find nearest non-missing neighbour for missing probes:
        nn <- findNN(pos=position,obs=obs)
     
        #Include probes with missing values in segments where their nearest neighbour probes are located
        new.res <- handleMissing(nn=nn,pos=position,obs=obs,pos.start=pos.start,pos.stop=pos.stop,seg.npos=seg.npos)  
        pos.start <- new.res$pos.start
        pos.stop <- new.res$pos.stop
        seg.npos <- new.res$seg.npos
       
        if(return.est){
          yhat[!obs] <- yhat[nn]
        }
			}
			
		}else{
		  warning(paste("pcf is not run for sample ",i," because all observations are missing. NA is returned.",sep=""),immediate.=TRUE,call.=FALSE)
		  seg.start <- 1
		  seg.stop <- nProbe
	    pos.start <- position[seg.start]
	    pos.stop <- position[seg.stop]
	    nSeg <- 1
	    seg.mean <- NA
	    seg.npos <- nProbe
	  }  
		
		
		#Round:
		seg.mean <- round(seg.mean,digits=digits)
		
    #Add results for this sample to results for other samples in data frame:
    seg <- data.frame(rep(sampleid[i],nSeg),pos.start,pos.stop,seg.npos,seg.mean,stringsAsFactors=FALSE)
		colnames(seg) <- seg.names
		segments <- rbind(segments,seg)

    if(return.est){
      #Rounding:
		  yhat <- round(yhat,digits=digits)
      pcf.est <- cbind(pcf.est,yhat)
    }
  }#endfor
  
  
    	
  if(verbose){
    cat(paste("pcf finished for sample ",i,sep=""),"\n")
  }

  
	if(return.est){
    pcf.est <- data.frame(position,pcf.est,stringsAsFactors=FALSE)
	  colnames(pcf.est) <- pcf.names
		return(list(estimates=pcf.est,segments=segments))
	}else{
	 return(segments)
	}
	
}#endfunction

