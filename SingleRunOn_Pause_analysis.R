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
                         dist=200, step=2, up_dist =20, file.pdf="heatmap.pdf", 
                         bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 1, use.log=FALSE, times=1){
  # AT here is the window of pause
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
  
  length_order  <- order(AT$long.pause - AT$short.paus, decreasing = T)
  AT <- AT[length_order  ,]
  hmat.pat.peak <-  hmat.pat.peak[length_order]
  hmat.mat.peak <- hmat.mat.peak[length_order]
  AT$V2 <- AT$V2 + AT$short.pause * step
  AT$V3 <- AT$V2 + AT$long.pause * step

  
  # make all beds the same length
  # plus strand chromEnd = chromStart + dist
  # minus strand chromStart = chromEnd - dist
  bed6 <- AT
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- AT[i,3] + up_dist
      bed6[i,2] <- AT[i,3] - dist}
    else
    {bed6[i,2] <- AT[i,2] - up_dist
    bed6[i,3] <- AT[i,2] + dist
    }
  }

  hmat.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, bed6[,c(1:6)],  step, times=times, use.log=use.log) 
  

  AT$start.steps <- up_dist%/%step+1
  AT$end.steps <- (AT$V3-AT$V2+up_dist)%/%step+1
  bin_number <- (up_dist + dist)/step
  
  # long pause
  hmat.high <- hmat.mat
  hmat.high[hmat.pat.peak >hmat.mat.peak, ] = hmat.pat[hmat.pat.peak >hmat.mat.peak, ] 
  
  # short pause
  hmat.low <- hmat.pat
  hmat.low[hmat.pat.peak >hmat.mat.peak, ] = hmat.mat[hmat.pat.peak >hmat.mat.peak, ]
  
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

    return (AT)
}



heatmap.Pause_allreads <-function(AT, file.plus.bw, file.minus.bw ,
                                  dist=200, step=2, up_dist =50, file.pdf="heatmap.pdf",
                                  bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 1, use.log=FALSE, times=1){
  AT <- heatmap.Pause(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat,
                  breaks=breaks, show.AT.line=FALSE,  navg = navg, times=times,
                  dist=dist, step=step, up_dist =up_dist , file.pdf=paste("allelicReads_heatmap_step",step,"_up",up_dist,"_dist",dist,".pdf",sep = ""))

  
  
  # make all beds the same length
  # plus strand chromEnd = chromStart + dist
  # minus strand chromStart = chromEnd - dist
  bed6 <- AT
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- AT[i,3] + up_dist
      bed6[i,2] <- AT[i,3] - dist}
    else
    {bed6[i,2] <- AT[i,2] - up_dist
    bed6[i,3] <- AT[i,2] + dist
    }
  }
  
  # get the sum of reads in each bin (size = step)
  hmat.pat <- read_read_mat_S (file.plus.bw, file.minus.bw, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <-hmat.pat 
  
  
  
  AT$start.steps <- up_dist%/%step+1
  AT$end.steps <- (AT$V3-AT$V2+up_dist)%/%step+1
  bin_number <- (up_dist + dist)/step
  
  
  hmat.high <- hmat.mat
  hmat.low <- hmat.mat
  
  
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

### SNP locations
# make SNPs location as bigwig files

heatmap.SNPsLocation.inPause <-function(AT, SNP.bw, file.plus.bw, file.minus.bw ,
                                        dist=200, step=2, up_dist =50, file.pdf="heatmap.pdf",
                                        bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 1, use.log=FALSE, times=1){
  AT <- heatmap.Pause(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat)
  
  # make all beds the same length
  # plus strand chromEnd = chromStart + dist
  # minus strand chromStart = chromEnd - dist
  bed6 <- AT
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- AT[i,3] + up_dist
      bed6[i,2] <- AT[i,3] - dist}
    else
    {bed6[i,2] <- AT[i,2] - up_dist
    bed6[i,3] <- AT[i,2] + dist
    }
  }
  
  # get the sum of reads in each bin (size = step)
  hmat.pat <- read_read_mat_SNPs (SNP.bw, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <-hmat.pat 
  
  
  
  AT$start.steps <- up_dist%/%step+1
  AT$end.steps <- (AT$V3-AT$V2+up_dist)%/%step+1
  bin_number <- (up_dist + dist)/step
  
  
  hmat.high <- hmat.mat
  hmat.low <- hmat.mat
  
  
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


t="KD"
pause_window <- read.table(paste("/Volumes/SPC_SD/KD_IGV/",t,"_dREG_5mat5pat_uniq_pValue_fdr0.1.bed.bed", sep =""), header = F)
AT <- pause_window
AT <- AT[AT$V1 != 'chrX',]
end=".rpm.bw"; times=10
end=".bw"; times=1
#HT.mat.map2ref.1bp_plus.bw
file.bw.plus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_plus",end, sep="")
file.bw.minus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".pat.map2ref.1bp_minus",end, sep="")
file.bw.plus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_plus",end, sep="")
file.bw.minus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,".mat.map2ref.1bp_minus",end, sep="")
file.plus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_plus",end, sep="")
file.minus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5N6_dedup_QC_end_minus",end, sep="")
SNP.bw <- "/Volumes/SPC_SD/KD_IGV/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"

breaks =seq(0, 10, 0.001)


if(0){
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
                       breaks=seq(0, 20, 0.001), navg = navg)
}

heatmap.Pause_allreads(AT, file.plus.bw, file.minus.bw ,
                       dist=100, step=2, up_dist =50, file.pdf="allreads_heatmap_step2_up50_dist100.pdf",
                       breaks=seq(0, 20, 0.001), navg = navg)


heatmap.Pause_allreads(AT, file.plus.bw, file.minus.bw ,
                       dist=200, step=10, up_dist =50, file.pdf="allreads_heatmap_step10_up50_dist200.pdf",
                       breaks=seq(0, 20, 0.001), navg = navg)



heatmap.SNPsLocation.inPause (AT, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=200, step=10, up_dist =50, file.pdf="SNPs_heatmap_step10_up50_dist200.pdf",
                              breaks=seq(0, 4, 0.1))
heatmap.SNPsLocation.inPause (AT, SNP.bw, file.plus.bw, file.minus.bw ,
                              dist=100, step=2, up_dist =50, file.pdf="SNPs_heatmap_step2_up50_dist100.pdf",
                              breaks=seq(0, 1, 0.1))
                       
