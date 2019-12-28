source("heatmaps.R");
source("draw-AllelelicTermination-heatmap.R");

getProfile <- function( bed1, hMarkFile, dist= 5000, step=25 )
{
	## Load mark.
	hMark <- load.bigWig(hMarkFile)
	mat1 <- NULL;

	if(NROW(bed1)!=0)
	{
		if(NROW(bed1)<100)
		{
			hCountMatrix <- bed.step.bpQuery.bigWig(hMark, center.bed(bed1[,c(1,2,3)], dist, dist), step=step, abs.value=TRUE)
			hmat1 <- matrix(unlist(hCountMatrix), nrow= NROW(bed1), byrow=TRUE);
			mat1 <- colMeans(hmat1);
		}
		else
		{
			hmat1 <- metaprofile.bigWig(center.bed(bed1[,c(1,2,3)], dist, dist), hMark, step=step)
			mat1 <- abs(hmat1$middle)
		}
	}

    unload.bigWig(hMark);

    return(mat1)
}


# library(matrixStats)

metaplot.AT <-function(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, 
                       name="", map5=TRUE, 
                       up_dist =50000, dist=50000, step=100,  file.pdf="heatmap.pdf", 
                       show.window=5000, 
                       ylim=c(0,1), 
                       xlab="", ylab= "proseq signal",
                       hl.bw.plus.pat=NULL,hl.bw.minus.pat=NULL, hl.bw.plus.mat=NULL ,hl.bw.minus.mat=NULL, 
                       use.log=FALSE){
  AT <- AT[,1:6]
  # make all beds the same length
  # plus strand chromEnd = chromStart + dist
  # minus strand chromStart = chromEnd - dist
  bed6 <- AT
  if (map5){
    for (i in 1:NROW(bed6)){
      if(bed6[i,6]=="-") {
        bed6[i,3] <- AT[i,3] + up_dist
        bed6[i,2] <- AT[i,3] - dist}
      else
      {bed6[i,2] <- AT[i,2] - up_dist
      bed6[i,3] <- AT[i,2] + dist
      }
    }
    AT$start.steps <- up_dist%/%step+1
    AT$end.steps <- (AT$V3-AT$V2+up_dist)%/%step+1
    bin_number <- (up_dist + dist)/step
  }else{
    for (i in 1:NROW(bed6)){
      if(bed6[i,6]=="-") {
        bed6[i,3] <- AT[i,2] + up_dist
        bed6[i,2] <- AT[i,2] - dist}
      else
      {bed6[i,2] <- AT[i,3] - up_dist
      bed6[i,3] <- AT[i,3] + dist
      }
    }
    AT$start.steps <- (up_dist - (AT$V3-AT$V2))%/%step + 1
    AT$end.steps <- up_dist%/%step + 1
    bin_number <- (up_dist + dist)/step
  }
  
  # get the sum of reads in each bin (size = step)
  hmat.pat <- read_read_mat2 (file.bw.plus.pat, file.bw.minus.pat, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <- read_read_mat2 (file.bw.plus.mat, file.bw.minus.mat, bed6[,c(1:6)],  step, times=times, use.log=use.log)     
  
  
  
  if(is.null(hl.bw.plus.pat)) {
    # dertemine High/Low allele based on the reads count within the AT regions (Before adjust by dist)
    hmat.pat.AT.rowSums <- NULL
    hmat.mat.AT.rowSums <- NULL
    for (i in 1:NROW(hmat.pat)){
      hmat.pat.AT.rowSums[i] <- sum(hmat.pat[i,][max(0,AT$start.steps[i]):min(AT$end.steps[i],bin_number)]) 
      hmat.mat.AT.rowSums[i] <- sum(hmat.mat[i,][max(0,AT$start.steps[i]):min(AT$end.steps[i],bin_number)]) 
    }
  } else {
    hl.pat <- read_read_mat2 (hl.bw.plus.pat, hl.bw.minus.pat, bed6[,c(1:6)], step, times=times)
    hl.mat  <- read_read_mat2 (hl.bw.plus.mat, hl.bw.minus.mat, bed6[,c(1:6)],  step, times=times)
    hmat.pat.AT.rowSums <- NULL
    hmat.mat.AT.rowSums <- NULL
    for (i in 1:NROW(hmat.pat)){
      hmat.pat.AT.rowSums[i] <- sum(hmat.pat[i,][max(0,AT$start.steps[i]):min(AT$end.steps[i],bin_number)]) 
      hmat.mat.AT.rowSums[i] <- sum(hmat.mat[i,][max(0,AT$start.steps[i]):min(AT$end.steps[i],bin_number)]) 
    }
  }
  hmat.high <- hmat.mat
  hmat.high[hmat.pat.AT.rowSums > hmat.mat.AT.rowSums, ] = hmat.pat[hmat.pat.AT.rowSums > hmat.mat.AT.rowSums, ] 
  
  hmat.low <- hmat.mat
  hmat.low[hmat.pat.AT.rowSums < hmat.mat.AT.rowSums, ] = hmat.pat[hmat.pat.AT.rowSums < hmat.mat.AT.rowSums, ] 
  
  if(NROW(bed6)<100){
    meta.hmat.high <- colMedians(hmat.high)
    meta.hmat.low <- colMedians(hmat.low)
  }else{
    meta.hmat.high.temp <- subsampled.quantiles.metaprofile(hmat.high)
    meta.hmat.high <- meta.hmat.high.temp$middle
    meta.hmat.low.temp <- subsampled.quantiles.metaprofile(hmat.low)
    meta.hmat.low <- meta.hmat.low.temp$middle
  }
  show.window <- min(show.window, (up_dist+dist))
  plot(seq(-1*(show.window)+step/2, show.window, step), meta.hmat.high[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ], col="red", 
       main = name,
       ylab= ylab,
       type = "l", 
       ylim = ylim,
       xlab=xlab
  )
  lines(seq(-1*(show.window)+step/2, show.window, step), meta.hmat.low[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ], col="blue")
  legend("topright", legend=c("high", "low"),
         col=c("red", "blue"), bty = "n", lty=1)
}


metaplot.AT.allReads <-function(AT, file.bw.plus,file.bw.minus, 
                                name="", map5=TRUE, add=FALSE, col="red", 
                                up_dist =50000, dist=50000, step=1000,  file.pdf="heatmap.pdf",
                                show.window=5000, 
                                ylim=c(0,20), lty=1,
                                xlab="", ylab= "proseq signal",
                                use.log=FALSE){
  AT <- AT[,1:6]
  # make all beds the same length
  # plus strand chromEnd = chromStart + dist
  # minus strand chromStart = chromEnd - dist
  bed6 <- AT
  if (map5){
    for (i in 1:NROW(bed6)){
      if(bed6[i,6]=="-") {
        bed6[i,3] <- AT[i,3] + up_dist
        bed6[i,2] <- AT[i,3] - dist}
      else
      {bed6[i,2] <- AT[i,2] - up_dist
      bed6[i,3] <- AT[i,2] + dist
      }
    }
    AT$start.steps <- up_dist%/%step+1
    AT$end.steps <- (AT$V3-AT$V2+up_dist)%/%step+1
    bin_number <- (up_dist + dist)/step
  }else{
    for (i in 1:NROW(bed6)){
      if(bed6[i,6]=="-") {
        bed6[i,3] <- AT[i,2] + up_dist
        bed6[i,2] <- AT[i,2] - dist}
      else
      {bed6[i,2] <- AT[i,3] - up_dist
      bed6[i,3] <- AT[i,3] + dist
      }
    }
    AT$start.steps <- (up_dist - (AT$V3-AT$V2))%/%step + 1
    AT$end.steps <- up_dist%/%step + 1
    bin_number <- (up_dist + dist)/step
  }
  
  # get the sum of reads in each bin (size = step)
  hmat.high <- read_read_mat2 (file.bw.plus, file.bw.minus, bed6[,c(1:6)], step, times=times, use.log=use.log)
  
  if(NROW(bed6)<100){
    meta.hmat.high <- colMedians(hmat.high)
  }else{
    meta.hmat.high.temp <- subsampled.quantiles.metaprofile(hmat.high)
    meta.hmat.high <- meta.hmat.high.temp$middle
  }
  
  show.window <- min(show.window, (up_dist+dist))
  if (add){
    lines(seq(-1*(show.window)+step/2, show.window, step), meta.hmat.high[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ], col=col, lty=lty)
    
  }else{
    plot(seq(-1*(show.window)+step/2, show.window, step), meta.hmat.high[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ], 
         col=col, 
         main = name,
         ylab= ylab,
         type = "l", 
         ylim = ylim,
         xlab=xlab,
         lty=lty
    )
  }
}