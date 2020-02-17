setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Kidney_and_SingleRunOn/Pause_manuscript")
source("/Users/shaopei/Box\ Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/heatmap/heatmaps.R")


draw_legend <- function(bk, hmcols, use.log=FALSE)
{
    par(mar=c(10,0,0,0), plt=c(0.15, 1.00, 0.4, 0.8 ));
    plot(NA,NA, type="n", xlim=c(1, NROW(bk)*0.85/0.7), ylim=c(0,1), xlab="", ylab="", bty="n", xaxt="n", yaxt= "n",bty="n" );
    for(i in 1:NROW(bk)) rect(i,0, i+1, 1, col=hmcols[i], border=hmcols[i]);
    if (use.log){
        axis(1, c(1, NROW(bk)/2, NROW(bk)), round(exp(c(bk[1], bk[round(NROW(bk)/2)], bk[NROW(bk)])),0), cex.axis=3, cex=3, tick=FALSE );
    } else {
        axis(1, c(1, NROW(bk)/2, NROW(bk)), round((c(bk[1], bk[round(NROW(bk)/2)], bk[NROW(bk)])),0), cex.axis=3, cex=3, tick=FALSE );
        }
}


read_read_mat_S <-function (file.plus.bw, file.minus.bw , bed6, step=2, navg = 20, times=1, use.log=FALSE)
{
    bw.plus  <- load.bigWig( file.plus.bw )
    bw.minus <- load.bigWig( file.minus.bw )

    hCountMatrix <- bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6[,c(1:6)] , step=step, abs.value=TRUE, op = "sum")
    hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(bed6[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
    if (!use.log){
        hmat <- times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) ;
    } else {
        hmat <- log(times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) + 1) ;
        }
    #avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
    
    unload.bigWig(bw.plus);
    unload.bigWig(bw.minus);
    
    return(hmat);    
}

read_read_mat_SNPs <-function (SNP.bw , bed6, step=2, navg = 20, times=1, use.log=FALSE)
{
  bw.plus  <- load.bigWig(SNP.bw )
  
  hCountMatrix <- bed.step.bpQuery.bigWig(bw.plus, bed6[,c(1:3)] , step=step, abs.value=TRUE, op = "sum")
  hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(bed6[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
  if (!use.log){
    hmat <- times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) ;
  } else {
    hmat <- log(times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) + 1) ;
  }
  #avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
  
  unload.bigWig(bw.plus);
  
  return(hmat);    
}


heatmap.Pause <-function(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat ,
                         dist=200, step=2, up_dist =20, file.pdf="heatmap.pdf", map5=TRUE, metaplot.pdf="metaplot.pdf",
                         heatmap=TRUE, metaplot=TRUE, metaplot.ylim=c(0,20),
                         bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=FALSE, navg = 1, use.log=FALSE, times=1){
  # AT here is the window of pause, dREG sites
  AT <- AT[,1:6] 
  
  # identify the location of max reads in mat and pat pause
  # determine the location of long pause and short pause
  hmat.pat.peak <- NULL
  hmat.mat.peak <- NULL
  for (i in 1:NROW(AT)){
    hmat.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, AT[i,], step, times=times, use.log=use.log)
    hmat.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, AT[i,], step, times=times, use.log=use.log) 
    
    hmat.pat.peak[i] <- which(hmat.pat == max(hmat.pat))
    hmat.mat.peak[i] <- which(hmat.mat == max(hmat.mat))
    AT$short.pause[i] = min(hmat.pat.peak[i],  hmat.mat.peak[i]) 
    AT$long.pause[i] = max(hmat.pat.peak[i],  hmat.mat.peak[i]) 
  }
  
  length_order  <- order(AT$long.pause - AT$short.pause, decreasing = T)
  AT$l = AT$long.pause - AT$short.pause
  AT <- AT[length_order  ,]
  hmat.pat.peak <-  hmat.pat.peak[length_order]
  hmat.mat.peak <- hmat.mat.peak[length_order]

  # new AT is the the wondow between short and long pause
  # + strand, V2=short, V3=long
  # - strand, V3=short, V2=long
  AT$old.V2 = AT$V2
  AT$old.V3 = AT$V3
  for (i in 1:NROW(AT)){
    if(AT[i,6]=="-") {
      AT$V2[i] <- AT$old.V3[i] - (AT$long.pause[i]) * step 
      AT$V3[i] <- AT$old.V3[i] - (AT$short.pause[i] -1)* step
    }
    else{
      AT$V3[i] <- AT$old.V2[i] + ((AT$long.pause[i]) * step)
      AT$V2[i] <- AT$old.V2[i] + (AT$short.pause[i] -1) * step
      
    }
  }
  

  # for heatmap, need to make all beds the same length
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
  
  hmat.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, bed6[,c(1:6)],  step, times=times, use.log=use.log) 
  
  
  
  # long pause
  hmat.high <- hmat.mat
  hmat.high[hmat.pat.peak >hmat.mat.peak, ] = hmat.pat[hmat.pat.peak >hmat.mat.peak, ] 
  
  # short pause
  hmat.low <- hmat.pat
  hmat.low[hmat.pat.peak >hmat.mat.peak, ] = hmat.mat[hmat.pat.peak >hmat.mat.peak, ]
  
  # metaplot
if(metaplot){  
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
  pdf(metaplot.pdf, width=6, height = 6)
  par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  par(mgp=c(3,1,0))
  par(cex.lab=2.2, cex.axis=2.2)
  plot(seq(-1*(show.window)+step/2, show.window, step), meta.hmat.high[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ], col="red", 
       main = "",
       ylab= "proseq signal",
       type = "l", 
       xlab="",
       ylim=metaplot.ylim,
       las=1
  )
  lines(seq(-1*(show.window)+step/2, show.window, step), meta.hmat.low[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ], col="blue")
  legend("topright", legend=c("late", "early"),
         col=c("red", "blue"), bty = "n", lty=1)
  dev.off()
}
  if (heatmap){
  
        if (navg > 1){
          hmat.low <- avgMat(hmat.low, navg = navg)
          hmat.high <- avgMat(hmat.high, navg = navg)
          newAT <- cbind(AT$start.steps, AT$end.steps)
          aveAT <- avgMat(newAT, navg = navg)
        }else {
          aveAT<- cbind(AT$start.steps, AT$end.steps)
        }

   #save.image(paste(file.pdf,".RData", sep=""))
   #load("data-hmat.RData")
   ## Write out a heatmap.
    if(is.null(breaks)) {
        bk.low <- seq(min(hmat.low), max(hmat.low), 0.01)
    } else {
        bk.low <- breaks
    }

    if(is.null(cols)) {
        hmcols.low <- (colorRampPalette(brewer.pal(9,"Reds"))(length(bk.low)-1))
    } else {
        hmcols.low <- colorRampPalette(cols)(length(bk.low)-1) # red
    }

    if(is.null(breaks)) {
        bk.high <- seq(min(hmat.high), max(hmat.high), 0.01)
    } else {
        bk.high <- breaks
    }

    if(is.null(cols)) {
        hmcols.high <- (colorRampPalette(brewer.pal(9,"Reds"))(length(bk.high)-1))
    } else {
        hmcols.high <- colorRampPalette(cols)(length(bk.high)-1) # red
    }


   # draw heatmap based on the score in the bins
   pdf(file.pdf, width=20, height = 20 )
   lay.heights <- c(0.85, 0.15);
   lay.widths  <- c(0.48, 0.04, 0.48 )
   layout(matrix(c(1, 5, 2, 3, 5, 4 ), nrow=2, byrow=T), widths=lay.widths, heights=lay.heights)

    ##part 1:
    gt <- pheatmap( hmat.low , cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.low, breaks = bk.low, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T, border_color = NA )
    par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
    ##grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
    gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
    # add a line to indicates the boundry of AT
    if (show.AT.line) {
        for (i in 1:NROW(hmat.high)){
          if (aveAT[i,2] <= bin_number){
          gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((floor(aveAT[i,2])-bl_wd),0) : floor(aveAT[i,2]) ] <- '#E69F00'; # "#88419D" ; #purple
          }
        }
    }
    grid.draw(gt$gtable)
    popViewport()

    ##part 2
    gt <- pheatmap( hmat.high, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.high, breaks = bk.high, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T, border_color = NA )
    par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
    plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
    ##grid.newpage()
    pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
    gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 3 );
    # add a line to indicates the boundry of AT
    if (show.AT.line) {
        for (i in 1:NROW(hmat.high)){
          if (aveAT[i,2] <= bin_number){
          gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((floor(aveAT[i,2])-bl_wd),0) : floor(aveAT[i,2]) ] <- '#E69F00'; #"#88419D"  ;
          }
        }
    }

    grid.draw(gt$gtable)
    popViewport()

    ##part 3:
    ## draw colorScale for peak track
    draw_legend(bk.low, hmcols.low, use.log =  use.log );
    
    ##part 4:
    ## draw colorScale for Heatmap
    draw_legend(bk.high, hmcols.high, use.log =  use.log  );
    dev.off()
  }

    return (AT)
}



heatmap.Pause_allreads <-function(AT, file.plus.bw, file.minus.bw ,
                                  dist=200, step=2, up_dist =50, file.pdf="heatmap.pdf", metaplot.pdf="metaplot.pdf", 
                                  allelic.heatmap.pdf="allelic.heatmap.pdf", allelic.metaplot.pdf="allelic.metaplot.pdf",
                                  bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 1, use.log=FALSE, times=1, 
                                  heatmap=TRUE, metaplot=TRUE, map5=TRUE, metaplot.ylim=c(0,30)){
  AT <- heatmap.Pause(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat,
                  breaks=breaks, show.AT.line=FALSE,  navg = navg, times=times, map5=map5,
                  dist=dist, step=step, up_dist =up_dist , metaplot.ylim = metaplot.ylim,
                  file.pdf= allelic.heatmap.pdf, metaplot.pdf= allelic.metaplot.pdf)

  
  
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
  hmat.pat <- read_read_mat_S (file.plus.bw, file.minus.bw, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <-hmat.pat 
  

  hmat.high <- hmat.mat
  hmat.low <- hmat.mat
  


  # metaplot
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
  pdf(metaplot.pdf, width=6, height = 6)
  par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  par(mgp=c(3,1,0))
  par(cex.lab=2.2, cex.axis=2.2)
  
  plot(seq(-1*(show.window)+step/2, show.window, step), meta.hmat.high[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ], col="purple", 
       main = "",
       ylab= "proseq signal",
       type = "l", 
       xlab="",
       ylim=metaplot.ylim,
       las=1
      
  )
#  lines(seq(-1*(show.window)+step/2, show.window, step), meta.hmat.low[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ], col="blue")
#  legend("topright", legend=c("late", "early"),
#         col=c("red", "blue"), bty = "n", lty=1)
  dev.off()
  #
  
  #heatmap
  if (navg > 1){
    hmat.low <- avgMat(hmat.low, navg = navg)
    hmat.high <- avgMat(hmat.high, navg = navg)
    newAT <- cbind(AT$start.steps, AT$end.steps)
    aveAT <- avgMat(newAT, navg = navg)
  }else {
    aveAT<- cbind(AT$start.steps, AT$end.steps)
  }
  
  if (heatmap){
  #save.image(paste(file.pdf,".RData", sep=""))
  #load("data-hmat.RData")
  ## Write out a heatmap.
  if(is.null(breaks)) {
    bk.low <- seq(min(hmat.low), max(hmat.low), 0.01)
  } else {
    bk.low <- breaks
  }
  
  if(is.null(cols)) {
    hmcols.low <- (colorRampPalette(brewer.pal(9,"Reds"))(length(bk.low)-1))
  } else {
    hmcols.low <- colorRampPalette(cols)(length(bk.low)-1) # red
  }
  
  if(is.null(breaks)) {
    bk.high <- seq(min(hmat.high), max(hmat.high), 0.01)
  } else {
    bk.high <- breaks
  }
  
  if(is.null(cols)) {
    hmcols.high <- (colorRampPalette(brewer.pal(9,"Reds"))(length(bk.high)-1))
  } else {
    hmcols.high <- colorRampPalette(cols)(length(bk.high)-1) # red
  }
  
  
  # draw heatmap based on the score in the bins
  pdf(file.pdf, width=20, height = 20 )
  lay.heights <- c(0.85, 0.15);
  lay.widths  <- c(0.48, 0.04, 0.48 )
  layout(matrix(c(1, 5, 2, 3, 5, 4 ), nrow=2, byrow=T), widths=lay.widths, heights=lay.heights)
  
  ##part 1:
  gt <- pheatmap( hmat.low , cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.low, breaks = bk.low, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T, border_color = NA )
  par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
  plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
  ##grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
  gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
  # add a line to indicates the boundry of AT
  # if (show.AT.line) {
  #     for (i in 1:NROW(hmat.high)){
  #       if (aveAT[i,2] <= bin_number){
  #       gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((floor(aveAT[i,2])-bl_wd),0) : floor(aveAT[i,2]) ] <- '#E69F00'; # "#88419D" ; #purple
  #       }
  #     }
  # }
  grid.draw(gt$gtable)
  popViewport()
  
  ##part 2
  gt <- pheatmap( hmat.high, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.high, breaks = bk.high, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T, border_color = NA )
  par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
  plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
  ##grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
  gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 3 );
  # add a line to indicates the boundry of AT
  if (show.AT.line) {
    for (i in 1:NROW(hmat.high)){
      if (aveAT[i,2] <= bin_number){
        gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((floor(aveAT[i,2])-bl_wd),0) : floor(aveAT[i,2]) ] <- '#E69F00'; #"#88419D"  ;
      }
    }
  }
  
  grid.draw(gt$gtable)
  popViewport()
  
  ##part 3:
  ## draw colorScale for peak track
  draw_legend(bk.low, hmcols.low, use.log =  use.log );
  
  ##part 4:
  ## draw colorScale for Heatmap
  draw_legend(bk.high, hmcols.high, use.log =  use.log  );
  dev.off()
  }
}

### SNP locations
# make SNPs location as bigwig files

heatmap.SNPsLocation.inPause <-function(AT, SNP.bw, file.plus.bw, file.minus.bw ,
                                        dist=200, step=2, up_dist =50, file.pdf="heatmap.pdf", metaplot.pdf="metaplot.pdf",
                                        bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 1, use.log=FALSE, times=1,
                                        map5 =TRUE,
                                        heatmap=TRUE, metaplot=TRUE, use.sum=FALSE) {
  AT <- heatmap.Pause(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, heatmap=FALSE,  metaplot=FALSE)
  
   # make all beds the same length
  # plus strand chromEnd = chromStart + dist
  # minus strand chromStart = chromEnd - dist
  bed6 <- AT

    if (map5){
    for (i in 1:NROW(bed6)){
      if(bed6[i,6]=="-") {
        bed6[i,3] <- AT[i,3] + up_dist
        bed6[i,2] <- AT[i,3]-1 - dist} #BED file is 0-based, left open coordinates. eg: [0,5) = {0,1,2,3,4} include 0, exclude 5
      else
      {bed6[i,2] <- AT[i,2] - up_dist
      bed6[i,3] <- AT[i,2]+1 + dist  #BED file is 0-based, left open coordinates. eg: [0,5) = {0,1,2,3,4} include 0, exclude 5
      }
    }
    AT$start.steps <- up_dist%/%step+1
    AT$end.steps <- (AT$V3-AT$V2+up_dist)%/%step+1
    bin_number <- (up_dist + dist)/step
  }else{
    for (i in 1:NROW(bed6)){
      if(bed6[i,6]=="-") {
        bed6[i,3] <- AT[i,2]+1 + up_dist  #BED file is 0-based, left open coordinates. eg: [0,5) = {0,1,2,3,4} include 0, exclude 5
        bed6[i,2] <- AT[i,2] - dist}
      else
      {bed6[i,2] <- AT[i,3]-1 - up_dist  #BED file is 0-based, left open coordinates. eg: [0,5) = {0,1,2,3,4} include 0, exclude 5
      bed6[i,3] <- AT[i,3] + dist
      }
    }
    AT$start.steps <- (up_dist - (AT$V3-AT$V2))%/%step + 1
    AT$end.steps <- up_dist%/%step + 1
    bin_number <- (up_dist + dist)/step
  }

  # get the sum of reads in each bin (size = step)
  hmat.pat <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <-hmat.pat 
  
  hmat.high <- hmat.mat
  hmat.low <- hmat.mat
  
  # metaplot
  if(0){
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
  a=meta.hmat.high
  }
  
  if (use.sum){
    a = colSums(hmat.high)
  }else{
    a = colMeans(hmat.high)
  }

    show.window <- min(show.window, (up_dist+dist))
  pdf(metaplot.pdf, width=6, height = 6)
  par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  par(mgp=c(3,1,0))
  par(cex.lab=2.2, cex.axis=2.2)
  
  plot(seq(-1*(show.window)+step/2, show.window, step), a[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ], col="red", 
       main = "",
       ylab= "SNPs counts median",
       type = "l", 
       xlab="",
       las=1
       #ylim = c(0,0.5)
  )

  dev.off()
  
  # heatmap
  if(heatmap){
  if (navg > 1){
    hmat.low <- avgMat(hmat.low, navg = navg)
    hmat.high <- avgMat(hmat.high, navg = navg)
    newAT <- cbind(AT$start.steps, AT$end.steps)
    aveAT <- avgMat(newAT, navg = navg)
  }else {
    aveAT<- cbind(AT$start.steps, AT$end.steps)
  }
  
  #save.image(paste(file.pdf,".RData", sep=""))
  #load("data-hmat.RData")
  ## Write out a heatmap.
  if(is.null(breaks)) {
    bk.low <- seq(min(hmat.low), max(hmat.low), 0.01)
  } else {
    bk.low <- breaks
  }
  
  if(is.null(cols)) {
    hmcols.low <- (colorRampPalette(brewer.pal(9,"Reds"))(length(bk.low)-1))
  } else {
    hmcols.low <- colorRampPalette(cols)(length(bk.low)-1) # red
  }
  
  if(is.null(breaks)) {
    bk.high <- seq(min(hmat.high), max(hmat.high), 0.01)
  } else {
    bk.high <- breaks
  }
  
  if(is.null(cols)) {
    hmcols.high <- (colorRampPalette(brewer.pal(9,"Reds"))(length(bk.high)-1))
  } else {
    hmcols.high <- colorRampPalette(cols)(length(bk.high)-1) # red
  }
  
  
  # draw heatmap based on the score in the bins
  pdf(file.pdf, width=20, height = 20 )
  lay.heights <- c(0.85, 0.15);
  lay.widths  <- c(0.48, 0.04, 0.48 )
  layout(matrix(c(1, 5, 2, 3, 5, 4 ), nrow=2, byrow=T), widths=lay.widths, heights=lay.heights)
  
  ##part 1:
  gt <- pheatmap( hmat.low , cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.low, breaks = bk.low, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T, border_color = NA )
  par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
  plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
  ##grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
  gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
  # add a line to indicates the boundry of AT
  # if (show.AT.line) {
  #     for (i in 1:NROW(hmat.high)){
  #       if (aveAT[i,2] <= bin_number){
  #       gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((floor(aveAT[i,2])-bl_wd),0) : floor(aveAT[i,2]) ] <- '#E69F00'; # "#88419D" ; #purple
  #       }
  #     }
  # }
  grid.draw(gt$gtable)
  popViewport()
  
  ##part 2
  gt <- pheatmap( hmat.high, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.high, breaks = bk.high, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T, border_color = NA )
  par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
  plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
  ##grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
  gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 3 );
  # add a line to indicates the boundry of AT
  if (show.AT.line) {
    for (i in 1:NROW(hmat.high)){
      if (aveAT[i,2] <= bin_number){
        gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((floor(aveAT[i,2])-bl_wd),0) : floor(aveAT[i,2]) ] <- '#E69F00'; #"#88419D"  ;
      }
    }
  }
  
  grid.draw(gt$gtable)
  popViewport()
  
  ##part 3:
  ## draw colorScale for peak track
  draw_legend(bk.low, hmcols.low, use.log =  use.log );
  
  ##part 4:
  ## draw colorScale for Heatmap
  draw_legend(bk.high, hmcols.high, use.log =  use.log  );
  dev.off()
  }
  return (list(seq(-1*(show.window)+step/2, show.window, step), a[((up_dist-show.window)%/% step +1) :((up_dist+show.window)%/% step) ]))
}


step=2
for(t in c("HT", "SK", "KD")){
show.window=50
#end=".rpm.bw"; times=10
  end=".bw"; times=1
  #HT.mat.map2ref.1bp_plus.bw
  file.bw.plus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_plus",end, sep="")
  file.bw.minus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_minus",end, sep="")
  file.bw.plus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_plus",end, sep="")
  file.bw.minus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_minus",end, sep="")
  file.plus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_plus",end, sep="")
  file.minus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_minus",end, sep="")
  SNP.bw <- "/Volumes/SPC_SD/KD_IGV/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"
  
  pause_window_0.1 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.1.bed.bed", sep =""), header = F)
  pause_window_0.1 <- pause_window_0.1[pause_window_0.1$V1 != 'chrX',]
  pause_window_0.9 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.9.bed", sep =""), header = F)
  pause_window_0.9 <- pause_window_0.9[pause_window_0.9$V1 != 'chrX',]
  #AT <- AT[AT$V1 != 'chrX',]
  t0=t
  t=paste(t0,"_fdr0.1",sep = "")
  a_0.1 = heatmap.SNPsLocation.inPause (pause_window_0.1, SNP.bw, file.plus.bw, file.minus.bw ,
                                dist=200, step=step, up_dist =100, 
                                file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                                metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                                breaks=seq(0, 1, 0.1), map5=TRUE,heatmap=FALSE)
  t=paste(t0,"_fdr0.1",sep = "")
  a_0.9 = heatmap.SNPsLocation.inPause (pause_window_0.9, SNP.bw, file.plus.bw, file.minus.bw ,
                                        dist=200, step=step, up_dist =100, 
                                        file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                                        metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                                        breaks=seq(0, 1, 0.1), map5=TRUE,heatmap=FALSE)
  
  #dev.off()
  pdf(paste("SNP_fdr0.1N0.9_compare_",t0,".pdf",sep = ""))
  par(mfcol=c(3,2))
  plot(a_0.1[[1]],a_0.1[[2]]/a_0.9[[2]], type="o", xlab="center short pause",ylab="fdr0.1 colMean(SNPs)/(fdr0.9 colMean(SNPs)", las=1, main="ratio")
  abline(h=1,col="green")
  abline(v=4,col="green")
  abline(v=-4,col="green")
  abline(v=-20,col="green")
  plot(a_0.9[[1]], a_0.9[[2]], col="blue", xlab="center short pause", ylab="SNPs mean", main=t0, type="o", ylim=c(0,max(a_0.9[[2]],a_0.1[[2]])), las=1)
  points(a_0.1[[1]], a_0.1[[2]],col="red", type="o")
  abline(v=4,col="green")
  abline(v=-4,col="green")
  abline(v=-20,col="green")
  legend("topright", legend=c("fdr<=0.1", "fdr>0.9"),
         col=c("red", "blue"), bty = "n", lty=1, pch=1)
  plot(a_0.1[[1]],a_0.1[[2]] - a_0.9[[2]], type="o", xlab="center short pause",ylab="(fdr0.1 colMean(SNPs)) - (fdr0.9 colMean(SNPs)", las=1, main="substract")
  abline(h=0, col="green")
  abline(v=4,col="green")
  abline(v=-4,col="green")
  abline(v=-20,col="green")
  #dev.off()
  
  t=paste(t0,"_fdr0.1",sep = "")
  a_0.1 = heatmap.SNPsLocation.inPause (pause_window_0.1, SNP.bw, file.plus.bw, file.minus.bw ,
                                        dist=200, step=step, up_dist =100, 
                                        file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                                        metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                                        breaks=seq(0, 1, 0.1), map5=FALSE,heatmap=FALSE)
  t=paste(t0,"_fdr0.1",sep = "")
  a_0.9 = heatmap.SNPsLocation.inPause (pause_window_0.9, SNP.bw, file.plus.bw, file.minus.bw ,
                                        dist=200, step=step, up_dist =100, 
                                        file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                                        metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                                        breaks=seq(0, 1, 0.1), map5=FALSE,heatmap=FALSE)
  
  #dev.off()
  #pdf(paste("SNP_fdr0.1N0.9_compare_longPause_",t0,".pdf",sep = ""))
  #par(mfrow=c(3,1))
  plot(a_0.1[[1]],a_0.1[[2]]/a_0.9[[2]], type="o", xlab="center long pause",ylab="fdr0.1 colMean(SNPs)/(fdr0.9 colMean(SNPs)", las=1, main="ratio")
  abline(h=1,col="green")
  abline(v=4,col="green")
  abline(v=-8,col="green")
  abline(v=-20,col="green")
  plot(a_0.9[[1]], a_0.9[[2]], col="blue", xlab="center long pause", ylab="SNPs mean", main=t0, type="o", ylim=c(0,max(a_0.9[[2]],a_0.1[[2]])), las=1)
  points(a_0.1[[1]], a_0.1[[2]],col="red", type="o")
  abline(v=4,col="green")
  abline(v=-8,col="green")
  abline(v=-20,col="green")
  legend("topright", legend=c("fdr<=0.1", "fdr>0.9"),
         col=c("red", "blue"), bty = "n", lty=1, pch=1)
  plot(a_0.1[[1]],a_0.1[[2]] - a_0.9[[2]], type="o", xlab="center long pause",ylab="(fdr0.1 colMean(SNPs)) - (fdr0.9 colMean(SNPs)", las=1, main="substract")
  abline(h=0, col="green")
  abline(v=4,col="green")
  abline(v=-8,col="green")
  abline(v=-20,col="green")
  dev.off()
}
  # fisher's exact test
  # short pause sites
step=2
p.value <- NULL
odds.ratio <- NULL
for(t in c("HT", "SK", "KD")){
  show.window=50
  #end=".rpm.bw"; times=10
  end=".bw"; times=1
  #HT.mat.map2ref.1bp_plus.bw
  file.bw.plus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_plus",end, sep="")
  file.bw.minus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_minus",end, sep="")
  file.bw.plus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_plus",end, sep="")
  file.bw.minus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_minus",end, sep="")
  file.plus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_plus",end, sep="")
  file.minus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_minus",end, sep="")
  SNP.bw <- "/Volumes/SPC_SD/KD_IGV/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"
  
  pause_window_0.1 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.1.bed.bed", sep =""), header = F)
  pause_window_0.1 <- pause_window_0.1[pause_window_0.1$V1 != 'chrX',]
  pause_window_0.9 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.9.bed", sep =""), header = F)
  pause_window_0.9 <- pause_window_0.9[pause_window_0.9$V1 != 'chrX',]
  t0=t
  t=paste(t0,"_fdr0.1",sep = "")
  a_0.1 = heatmap.SNPsLocation.inPause (pause_window_0.1, SNP.bw, file.plus.bw, file.minus.bw ,
                                        dist=200, step=step, up_dist =100, 
                                        file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                                        metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                                        breaks=seq(0, 1, 0.1), map5=TRUE,heatmap=FALSE, use.sum = TRUE)
  t=paste(t0,"_fdr0.1",sep = "")
  a_0.9 = heatmap.SNPsLocation.inPause (pause_window_0.9, SNP.bw, file.plus.bw, file.minus.bw ,
                                        dist=200, step=step, up_dist =100, 
                                        file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                                        metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                                        breaks=seq(0, 1, 0.1), map5=TRUE,heatmap=FALSE, use.sum = TRUE)
  inside = which(a_0.1[[1]] > -4 & a_0.1[[1]] < 4 )
  outside = which(a_0.1[[1]] > -20 & a_0.1[[1]] < -4)
  
  testor = rbind(c(sum(a_0.1[[2]][inside]),sum(a_0.1[[2]][outside])),
              + c(sum(a_0.9[[2]][inside]),sum(a_0.9[[2]][outside])) ); testor
  f = fisher.test(testor); 
  p.value= c(p.value, f$p.value)
  odds.ratio = c(odds.ratio, f$estimate)
}
options(scipen = 1)
options(digits = 2)
p.value
odds.ratio
p.adjust(p.value, method="bonferroni")

# long pause sites
step=2
p.value <- NULL
odds.ratio <- NULL
for(t in c("HT", "SK", "KD")){
  show.window=50
  #end=".rpm.bw"; times=10
  end=".bw"; times=1
  #HT.mat.map2ref.1bp_plus.bw
  file.bw.plus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_plus",end, sep="")
  file.bw.minus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_minus",end, sep="")
  file.bw.plus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_plus",end, sep="")
  file.bw.minus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_minus",end, sep="")
  file.plus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_plus",end, sep="")
  file.minus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_minus",end, sep="")
  SNP.bw <- "/Volumes/SPC_SD/KD_IGV/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"
  
  pause_window_0.1 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.1.bed.bed", sep =""), header = F)
  pause_window_0.1 <- pause_window_0.1[pause_window_0.1$V1 != 'chrX',]
  pause_window_0.9 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.9.bed", sep =""), header = F)
  pause_window_0.9 <- pause_window_0.9[pause_window_0.9$V1 != 'chrX',]
  t0=t
  t=paste(t0,"_fdr0.1",sep = "")
  a_0.1 = heatmap.SNPsLocation.inPause (pause_window_0.1, SNP.bw, file.plus.bw, file.minus.bw ,
                                        dist=200, step=step, up_dist =100, 
                                        file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                                        metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                                        breaks=seq(0, 1, 0.1), map5=FALSE,heatmap=FALSE, use.sum = TRUE)
  t=paste(t0,"_fdr0.1",sep = "")
  a_0.9 = heatmap.SNPsLocation.inPause (pause_window_0.9, SNP.bw, file.plus.bw, file.minus.bw ,
                                        dist=200, step=step, up_dist =100, 
                                        file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                                        metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                                        breaks=seq(0, 1, 0.1), map5=FALSE,heatmap=FALSE, use.sum = TRUE)
  inside = which(a_0.1[[1]] > -8 & a_0.1[[1]] < 4 )
  outside = which(a_0.1[[1]] > -20 & a_0.1[[1]] < -8)
  
  testor = rbind(c(sum(a_0.1[[2]][inside]),sum(a_0.1[[2]][outside])),
                 + c(sum(a_0.9[[2]][inside]),sum(a_0.9[[2]][outside])) ); testor
  f = fisher.test(testor); f
  p.value= c(p.value, f$p.value)
  odds.ratio = c(odds.ratio, f$estimate)
}
p.value
odds.ratio
p.adjust(p.value, method="bonferroni")

# use heatmap to examine the proseq signal

t="HT"
show.window=100
for(t in c("HT", "SK", "KD")){
pause_window <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.1.bed.bed", sep =""), header = F)
#pause_window <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.9.bed", sep =""), header = F)
pause_window <- pause_window[pause_window$V1 != 'chrX',]

pause_window_0.1 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.1.bed.bed", sep =""), header = F)
pause_window_0.1 <- pause_window_0.1[pause_window_0.1$V1 != 'chrX',]
pause_window_0.9 <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.9.bed", sep =""), header = F)
pause_window_0.9 <- pause_window_0.9[pause_window_0.9$V1 != 'chrX',]

#AT <- AT[AT$V1 != 'chrX',]
end=".rpm.bw"; times=10
#end=".bw"; times=1
#HT.mat.map2ref.1bp_plus.bw
# allelic reads
file.bw.plus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_plus",end, sep="")
file.bw.minus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_minus",end, sep="")
file.bw.plus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_plus",end, sep="")
file.bw.minus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_minus",end, sep="")
# all reads
file.plus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_plus",end, sep="")
file.minus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_minus",end, sep="")
SNP.bw <- "/Volumes/SPC_SD/KD_IGV/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"

t0=t
t=paste(t0,"_fdr0.1",sep = "")
# ChRO-seq signal
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100, 
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160))

heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100, 
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map3.pdf",sep = ""),
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map3.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map3.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map3.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, heatmap = FALSE, map5=FALSE, metaplot.ylim = c(0,160))

# SNPs
heatmap.SNPsLocation.inPause (pause_window_0.1, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=200, step=2, up_dist =100, 
                              file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                              metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                              breaks=seq(0, 1, 0.1), map5=TRUE)


t=paste(t0,"_fdr0.9",sep = "")
heatmap.Pause_allreads(pause_window_0.9, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100, 
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map5.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map5.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE, metaplot.ylim = c(0,160))

heatmap.Pause_allreads(pause_window_0.9, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =100, 
                       file.pdf=paste(t,"_allreads_heatmap_step2_up100_dist200_map3.pdf",sep = ""),
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up100_dist200_map3.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up100_dist200_map3.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up100_dist200_map3.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, heatmap = FALSE, map5=FALSE, metaplot.ylim = c(0,160))

heatmap.SNPsLocation.inPause (pause_window_0.9, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=200, step=2, up_dist =100, 
                              file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                              metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                              breaks=seq(0, 1, 0.1), map5=TRUE)
}


# old scripts, not using
if(0){

t=paste(t0,"_fdr0.1",sep = "")
heatmap.SNPsLocation.inPause (pause_window_0.1, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=200, step=2, up_dist =100, 
                              file.pdf=paste(t,"_SNPs_heatmap_step2_up100_dist200_map5.pdf",sep=""),
                              metaplot.pdf=paste(t,"_SNPs_sum_step2_up100_dist200_map5.pdf",sep=""),
                              breaks=seq(0, 1, 0.1), map5=TRUE)

  # centher  at long pause (late)
heatmap.SNPsLocation.inPause (pause_window, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=100, step=2, up_dist =200, 
                              file.pdf=paste(t,"_SNPs_heatmap_step2_up200_dist100_map3.pdf",sep=""),
                              metaplot.pdf=paste(t,"_SNPs_sum_step2_up200_dist100_map3.pdf",sep=""),
                              breaks=seq(0, 1, 0.1), map5=FALSE)


t=paste(t0,"_fdr0.1",sep = "")
heatmap.Pause_allreads(pause_window_0.1, file.plus.bw, file.minus.bw ,
                       dist=100, step=2, up_dist =200, 
                       file.pdf=paste(t,"_allreads_heatmap_step2_up200_dist100_map3.pdf",sep = ""),
                       metaplot.pdf = paste(t,"_allreads_metaplot_step2_up200_dist100_map3.pdf",sep = ""),
                       allelic.heatmap.pdf = paste(t,"_allelicReads_heatmap_step2_up200_dist100_map3.pdf",sep = ""),
                       allelic.metaplot.pdf = paste(t,"_allelicReads__metaplot_step2_up200_dist100_map3.pdf",sep = ""),
                       breaks=seq(0, 20, 0.001), times=times, map5=FALSE)

  heatmap.Pause(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat,
                breaks=breaks)
heatmap.Pause(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat,
              breaks=breaks, show.AT.line=FALSE,
              dist=200, step=2, up_dist =50, file.pdf="allelicReads_heatmap_step2_up50_dist200.pdf")


#heatmap.Pause(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat,
#              breaks=breaks,
#              dist=1000, step=2, up_dist =500, file.pdf="heatmap_step2_up500_dist1000.pdf")




heatmap.Pause_allreads(AT, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =50, file.pdf="allreads_heatmap_step2_up50_dist200.pdf",
                       breaks=seq(0, 20, 0.001), times=times)

heatmap.Pause_allreads(AT, file.plus.bw, file.minus.bw ,
                       dist=100, step=2, up_dist =50, file.pdf="allreads_heatmap_step2_up50_dist100.pdf",
                       breaks=seq(0, 20, 0.001), times=times, map5=FALSE)


heatmap.Pause_allreads(AT, file.plus.bw, file.minus.bw ,
                       dist=200, step=10, up_dist =50, file.pdf="allreads_heatmap_step10_up50_dist200.pdf",
                       breaks=seq(0, 20, 0.001), times=times)

heatmap.Pause_allreads(AT, file.plus.bw, file.minus.bw ,
                       dist=50, step=10, up_dist =200, file.pdf="allreads_heatmap_step10_up50_dist200_map3.pdf",
                       breaks=seq(0, 20, 0.001), times=times, map5=FALSE)

heatmap.SNPsLocation.inPause (AT, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=200, step=10, up_dist =50, file.pdf="SNPs_heatmap_step10_up50_dist200.pdf",
                              breaks=seq(0, 4, 0.1), map5=TRUE)
heatmap.SNPsLocation.inPause (AT, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=100, step=2, up_dist =200, file.pdf="SNPs_heatmap_step2_up100_dist100.pdf",
                              breaks=seq(0, 1, 0.1), map5=TRUE)
heatmap.SNPsLocation.inPause (AT, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=100, step=2, up_dist =200, file.pdf="SNPs_heatmap_step2_up100_dist100_3prime.pdf",
                              breaks=seq(0, 1, 0.1), map5=FALSE)
heatmap.Pause_allreads(AT, file.plus.bw, file.minus.bw ,
                       dist=200, step=2, up_dist =200, file.pdf="allreads_heatmap_step2_up200_dist200_map5.pdf",
                       breaks=seq(0, 20, 0.001), times=times, map5=TRUE)
heatmap.SNPsLocation.inPause (AT, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=200, step=2, up_dist =200, file.pdf="SNPs_heatmap_step2_up200_dist200.pdf",
                              breaks=seq(0, 1, 0.1), map5=TRUE)
heatmap.SNPsLocation.inPause (AT, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=200, step=2, up_dist =200, file.pdf="SNPs_heatmap_step2_up200_dist200_map3.pdf",
                              breaks=seq(0, 1, 0.1), map5=FALSE)
}

