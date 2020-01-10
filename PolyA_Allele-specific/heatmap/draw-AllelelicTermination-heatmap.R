#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/heatmap")
source("heatmaps.R");
source("hist.param.R");
source("hist.svm.com.R");

### Load the package or install if not present
#if (!require("RColorBrewer")) {
#install.packages("RColorBrewer")
#library(RColorBrewer)
#}
#library(pheatmap)

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


read_read_mat0 <-function (file.bw.org, dregX, step, navg = 20, times=1)
{
    hmark <- load.bigWig( file.bw.org )

    hCountMatrix <- bed.step.bpQuery.bigWig( hmark, dregX[,c(1,2,3)], step=step, abs.value=TRUE)
    hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(dregX[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
    hmat <- log( times * matrix(unlist(hCountMatrix), nrow= NROW(dregX), byrow=TRUE)+1) ;

    avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
    
    unload.bigWig(hmark);
    
    return(avgMat);    
}

read_read_mat1 <-function (file.bw.org, dregX, step, navg = 1, times=1)
{
    hmark <- load.bigWig( file.bw.org )

    hCountMatrix <- bed.step.bpQuery.bigWig( hmark, dregX[,c(1,2,3)] , step=step, abs.value=TRUE)
    hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(dregX[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
    hmat <- times * matrix(unlist(hCountMatrix), nrow= NROW(dregX), byrow=TRUE) ;

    #avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
    
    unload.bigWig(hmark);
    
    return(hmat);    
}

read_read_mat2 <-function (file.plus.bw, file.minus.bw , bed6, step, navg = 20, times=1, use.log=FALSE)
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

read_read_mat2_log <-function (file.plus.bw, file.minus.bw , bed6, step, navg = 20, times=1)
{
    read_read_mat2(file.plus.bw, file.minus.bw , bed6, step, navg = navg, times=times, use.log=TRUE)
}

heatmap.gene<-function( df.bed.strand, file.plus.bw, file.minus.bw, file.bw.org, file.peak.org, file.bw.pred, file.peak.pred, file.pdf, subs = NULL, breaks = NULL, cols = NULL, step = 25) {
  if(is.na(file.bw.org) || is.na(file.bw.pred))
    return;
  
  dregX <- df.bed.strand 
  dregX <- dregX[grep("_|chrY|chrX", dregX[,1], invert=TRUE),]
  
  hmat.pred <- read_read_mat1 (file.bw.pred, dregX[,c(1:6)], step, times=1)
  hmat.org  <- read_read_mat1 (file.bw.org, dregX[,c(1:6)],  step, times=1)
  
  ## Write out a heatmap.
  if(is.null(breaks)) {
    bk.org <- seq(min(hmat.org), max(hmat.org), 1)
  } else {
    bk.org <- breaks
  }
  
  if(is.null(cols)) {
    hmcols.org <- rev(colorRampPalette(brewer.pal(9,"Reds"))(length(bk.org)-1))
  } else {
    hmcols.org <- colorRampPalette(cols)(length(bk.org)-1) # red
  }
  
  if(is.null(breaks)) {
    bk.pred <- seq(min(hmat.pred), max(hmat.pred), 1)
  } else {
    bk.pred <- breaks
  }
  
  if(is.null(cols)) {
    hmcols.pred <- rev(colorRampPalette(brewer.pal(9,"Reds"))(length(bk.pred)-1))
  } else {
    hmcols.pred <- colorRampPalette(cols)(length(bk.pred)-1) # red
  }
  
  ## start new pdf file
  pdf(file.pdf, width=20, height = 20 )
  
  
  lay.heights <- c(0.85, 0.15);
  lay.widths  <- c(0.48, 0.04, 0.48 )
  layout(matrix(c(1, 5, 2, 3, 5, 4 ), nrow=2, byrow=T), widths=lay.widths, heights=lay.heights)
  
  ##part2:
  gt <- pheatmap( hmat.pred, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.pred, breaks = bk.pred, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
  par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
  plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
  ##grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
  gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 3 );
  grid.draw(gt$gtable)
  popViewport()
  
  ##part 1:
  gt <- pheatmap( hmat.org, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols.org, breaks = bk.org, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
  par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
  plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
  ##grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
  gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
  grid.draw(gt$gtable)
  popViewport()
  
  ##part 3:
  ## draw colorScale for peak track
  draw_legend(bk.org, hmcols.org );
  
  ##part 4:
  ## draw colorScale for Heatmap
  draw_legend(bk.pred, hmcols.pred );
  
  
  dev.off(); 
  
  invisible()
  
} 

heatmap.AT0 <-function(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat ,dist, step, high.low.by=NULL, sort.by.HL.start = FALSE, file.pdf="heatmap.pdf", bl_wd=1, breaks=seq(0,10,1)){
  # high.low.by is bed regions
  
  AT <- AT[,1:6]
  if (sort.by.HL.start) {
    length_order  <- order((high.low.by$V2 - AT$V2), decreasing = T)
  } else {
    length_order  <- order(AT$V3 - AT$V2, decreasing = T)
  }
  AT <- AT[length_order  ,]
  # make all beds the same length
  # plus strand chromEnd = chromStart + dist
  # minus strand chromStart = chromEnd - dist
  bed6 <- AT
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,2] <- bed6[i,3]- dist}
    else
    {bed6[i,3] <- bed6[i,2] +dist
    }
  }
  
  # get the sum of reads in each bin (size = step)
  hmat.pat <- read_read_mat2 (file.bw.plus.pat, file.bw.minus.pat, bed6[,c(1:6)], step, times=1)
  hmat.mat  <- read_read_mat2 (file.bw.plus.mat, file.bw.minus.mat, bed6[,c(1:6)],  step, times=1)
  AT$steps <- (AT$V3-AT$V2)%/%step+1
  bin_number <- dist/step
  
  if(is.null(high.low.by)) {
    # dertemine High/Low allele based on the reads count within the AT regions (Before adjust by dist)
    hmat.pat.AT.rowSums <- NULL
    hmat.mat.AT.rowSums <- NULL
    for (i in 1:NROW(hmat.pat)){
      hmat.pat.AT.rowSums[i] <- sum(hmat.pat[i,][1:min(AT$steps[i],bin_number)]) 
      hmat.mat.AT.rowSums[i] <- sum(hmat.mat[i,][1:min(AT$steps[i],bin_number)]) 
    }
  } else {
    HL <- high.low.by[length_order  ,]
    HL$start <- (HL$V2-AT$V2)%/%step
    HL$end <- (HL$V3-AT$V2)%/%step+1
    hmat.pat.AT.rowSums <- NULL
    hmat.mat.AT.rowSums <- NULL
    for (i in 1:NROW(hmat.pat)){
      hmat.pat.AT.rowSums[i] <- sum(hmat.pat[i,][max(1,min(HL$start[i],bin_number)) : min(HL$end[i],bin_number)]) 
      hmat.mat.AT.rowSums[i] <- sum(hmat.mat[i,][max(1,min(HL$start[i],bin_number))  : min(HL$end[i],bin_number)]) 
    }
  }
  hmat.high <- hmat.mat
  hmat.high[hmat.pat.AT.rowSums > hmat.mat.AT.rowSums, ] = hmat.pat[hmat.pat.AT.rowSums > hmat.mat.AT.rowSums, ] 
  
  hmat.low <- hmat.mat
  hmat.low[hmat.pat.AT.rowSums < hmat.mat.AT.rowSums, ] = hmat.pat[hmat.pat.AT.rowSums < hmat.mat.AT.rowSums, ] 
  
  
  #save.image("data-hmat.RData")
  #load("data-hmat.RData")
  hmcols <- rev(colorRampPalette(brewer.pal(9,"RdBu"))(length(breaks)-1))
  # draw heatmap based on the score in the bins
  pdf(file.pdf, width=20, height = 20 )
  lay.heights <- c(0.85, 0.15);
  lay.widths  <- c(0.48, 0.04, 0.48 )
  layout(matrix(c(1, 5, 2, 3, 5, 4 ), nrow=2, byrow=T), widths=lay.widths, heights=lay.heights)
  
  ##part 1:
  gt <- pheatmap( hmat.low , cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = breaks, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
  par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
  plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
  ##grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
  gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 1)
  # add a black line to indicates the boundry of AT
  if(is.null(high.low.by)) {
    for (i in 1:NROW(hmat.high)){
      if (AT$steps[i] <= bin_number){
        gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((AT$steps[i]-bl_wd),0) : AT$steps[i] ] <- "#000000";
      }
    }
  } else {
    for (i in 1:NROW(hmat.high)){
      if (AT$steps[i] <= bin_number){
        gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((AT$steps[i]-bl_wd),0) : AT$steps[i] ] <- "#000000";
        gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((min(HL$start[i],bin_number)-bl_wd),0) : min(HL$start[i],bin_number) ] <- "#000000";
      }
    }
  }
  grid.draw(gt$gtable)
  popViewport()
  
  ##part 2
  gt <- pheatmap( hmat.high, cluster_rows = FALSE, cluster_cols = FALSE, col= hmcols, breaks = breaks, legend=FALSE, show_rownames=FALSE, show_colnames=FALSE, silent=T )
  par(mar=c(0,0,0,0), plt=c(0.2, 0.8,0.2, 0.8 ));
  plot(NA,NA, type="n", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1), xaxs="i", yaxs="i", xaxt = "n", yaxt = "n", bty="n");
  ##grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 3, widths=lay.widths, heights=lay.heights) ))
  gt$gtable$vp <- viewport(layout.pos.row = 1, layout.pos.col = 3 );
  # add a black line to indicates the boundry of AT
  if(is.null(high.low.by)) {
    for (i in 1:NROW(hmat.high)){
      if (AT$steps[i] <= bin_number){
        gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((AT$steps[i]-bl_wd),0) : AT$steps[i] ] <- "#000000";
      }
    }
  } else {
    for (i in 1:NROW(hmat.high)){
      if (AT$steps[i] <= bin_number){
        gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((AT$steps[i]-bl_wd),0) : AT$steps[i] ] <- "#000000";
        gt$gtable$grobs[[1]]$children[[1]]$gp$fill[i, max((min(HL$start[i],bin_number)-bl_wd),0) : min(HL$start[i],bin_number) ] <- "#000000";
      }
    }
  }
  grid.draw(gt$gtable)
  popViewport()
  
  ##part 3:
  ## draw colorScale for peak track
  draw_legend(breaks, hmcols );
  
  ##part 4:
  ## draw colorScale for Heatmap
  draw_legend(breaks, hmcols );
  dev.off()
  
}

avgMat <-function (hmat ,navg = 20){
    avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
}

heatmap.AT3 <-function(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat ,dist, step, up_dist =20000, file.pdf="heatmap.pdf", hl.bw.plus.pat=NULL,hl.bw.minus.pat=NULL, hl.bw.plus.mat=NULL ,hl.bw.minus.mat=NULL, bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 10, use.log=FALSE, times=1){
  AT <- AT[,1:6]
  length_order  <- order(AT$V3 - AT$V2, decreasing = T)
  AT <- AT[length_order  ,]
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
  hmat.pat <- read_read_mat2 (file.bw.plus.pat, file.bw.minus.pat, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <- read_read_mat2 (file.bw.plus.mat, file.bw.minus.mat, bed6[,c(1:6)],  step, times=times, use.log=use.log)     
  
  
  
  AT$start.steps <- up_dist%/%step+1
  AT$end.steps <- (AT$V3-AT$V2+up_dist)%/%step+1
  bin_number <- (up_dist + dist)/step
    
    if(is.null(hl.bw.plus.pat)) {
        # dertemine High/Low allele based on the reads count within the AT regions (Before adjust by dist)
        hmat.pat.AT.rowSums <- NULL
        hmat.mat.AT.rowSums <- NULL
        for (i in 1:NROW(hmat.pat)){
          hmat.pat.AT.rowSums[i] <- sum(hmat.pat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
          hmat.mat.AT.rowSums[i] <- sum(hmat.mat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
          }
    } else {
        hl.pat <- read_read_mat2 (hl.bw.plus.pat, hl.bw.minus.pat, bed6[,c(1:6)], step, times=times)
        hl.mat  <- read_read_mat2 (hl.bw.plus.mat, hl.bw.minus.mat, bed6[,c(1:6)],  step, times=times)
        hmat.pat.AT.rowSums <- NULL
        hmat.mat.AT.rowSums <- NULL
        for (i in 1:NROW(hmat.pat)){
          hmat.pat.AT.rowSums[i] <- sum(hl.pat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
          hmat.mat.AT.rowSums[i] <- sum(hl.mat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
          }
    }
        hmat.high <- hmat.mat
        hmat.high[hmat.pat.AT.rowSums > hmat.mat.AT.rowSums, ] = hmat.pat[hmat.pat.AT.rowSums > hmat.mat.AT.rowSums, ] 

        hmat.low <- hmat.mat
        hmat.low[hmat.pat.AT.rowSums < hmat.mat.AT.rowSums, ] = hmat.pat[hmat.pat.AT.rowSums < hmat.mat.AT.rowSums, ] 

        if (navg > 1){
          hmat.low <- avgMat(hmat.low, navg = navg)
          hmat.high <- avgMat(hmat.high, navg = navg)
          newAT <- cbind(AT$start.steps, AT$end.steps)
          aveAT <- avgMat(newAT, navg = navg)
        }
        else {
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

heatmap.AT4_allreads <-function(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat ,dist, step, up_dist =20000, file.pdf="heatmap.pdf", hl.bw.plus.pat=NULL,hl.bw.minus.pat=NULL, hl.bw.plus.mat=NULL ,hl.bw.minus.mat=NULL, bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 10, use.log=FALSE, times=1){
  AT <- AT[,1:6]
  length_order  <- order(AT$V3 - AT$V2, decreasing = T)
  AT <- AT[length_order  ,]
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
  hmat.pat <- read_read_mat2 (file.bw.plus.pat, file.bw.minus.pat, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <- read_read_mat2 (file.bw.plus.mat, file.bw.minus.mat, bed6[,c(1:6)],  step, times=times, use.log=use.log)     
  
  
  
  AT$start.steps <- up_dist%/%step+1
  AT$end.steps <- (AT$V3-AT$V2+up_dist)%/%step+1
  bin_number <- (up_dist + dist)/step
    
    if(is.null(hl.bw.plus.pat)) {
        # dertemine High/Low allele based on the reads count within the AT regions (Before adjust by dist)
        hmat.pat.AT.rowSums <- NULL
        hmat.mat.AT.rowSums <- NULL
        for (i in 1:NROW(hmat.pat)){
          hmat.pat.AT.rowSums[i] <- sum(hmat.pat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
          hmat.mat.AT.rowSums[i] <- sum(hmat.mat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
          }
    } else {
        hl.pat <- read_read_mat2 (hl.bw.plus.pat, hl.bw.minus.pat, bed6[,c(1:6)], step, times=times)
        hl.mat  <- read_read_mat2 (hl.bw.plus.mat, hl.bw.minus.mat, bed6[,c(1:6)],  step, times=times)
        hmat.pat.AT.rowSums <- NULL
        hmat.mat.AT.rowSums <- NULL
        for (i in 1:NROW(hmat.pat)){
          hmat.pat.AT.rowSums[i] <- sum(hl.pat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
          hmat.mat.AT.rowSums[i] <- sum(hl.mat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
          }
    }
        hmat.high <- hmat.mat
        hmat.high[hmat.pat.AT.rowSums > hmat.mat.AT.rowSums, ] = hmat.pat[hmat.pat.AT.rowSums > hmat.mat.AT.rowSums, ] 

        hmat.low <- hmat.mat
        hmat.low[hmat.pat.AT.rowSums < hmat.mat.AT.rowSums, ] = hmat.pat[hmat.pat.AT.rowSums < hmat.mat.AT.rowSums, ] 

        if (navg > 1){
          hmat.low <- avgMat(hmat.low, navg = navg)
          hmat.high <- avgMat(hmat.high, navg = navg)
          newAT <- cbind(AT$start.steps, AT$end.steps)
          aveAT <- avgMat(newAT, navg = navg)
        }
        else {
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


between_organs_AT_bw_heatmap_2 <- function (at, t){
    # use Alleleic Termination from one organ and see the proseq reads abundance from another organ
AT_intersect <- read.table(paste("../",at,"_AT_3tunitIntersectNativeHMM_intersectRegion.bed", sep=""), header = F)
file.bw.plus.pat <- paste("../",t,"_MB6_all_R1.pat_1bp_plus.bw", sep="")
file.bw.minus.pat <- paste("../",t,"_MB6_all_R1.pat_1bp_minus.bw", sep="")
file.bw.plus.mat <- paste("../",t,"_MB6_all_R1.mat_1bp_plus.bw", sep="")
file.bw.minus.mat <- paste("../",t,"_MB6_all_R1.mat_1bp_minus.bw", sep="")
heatmap.AT(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=50000, step=500,  file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT50Kb_step500.pdf",sep=""), bl_wd=1, breaks=seq(0,10,1))
heatmap.AT(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=100000, step=1000, file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT100Kb_step1K.pdf",sep="") , bl_wd=1, breaks=seq(0,10,1))
#heatmap.AT(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=1000000, step=1000, file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT1Mb_step1K.pdf",sep="") , bl_wd=1, breaks=seq(0,10,1))
heatmap.AT(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=200000, step=1000,file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT200Kb_step1K.pdf",sep="") , bl_wd=1, breaks=seq(0,10,1))
}

# between_organs_AT_bw_heatmap_2("BN","BN")
# between_organs_AT_bw_heatmap_2("LV","LV")
# between_organs_AT_bw_heatmap_2("LV","BN")
# between_organs_AT_bw_heatmap_2("BN","LV")

between_organs_AT_bw_heatmap_3 <- function (at, t, navg=5){
    # use Alleleic Termination from one organ and see the proseq reads abundance from another organ
AT_intersect <- read.table(paste("../",at,"_AT_3tunitIntersectNativeHMM_intersectRegion.bed", sep=""), header = F)
file.bw.plus.pat <- paste("../",t,"_MB6_all_R1.pat_1bp_plus.bw", sep="")
file.bw.minus.pat <- paste("../",t,"_MB6_all_R1.pat_1bp_minus.bw", sep="")
file.bw.plus.mat <- paste("../",t,"_MB6_all_R1.mat_1bp_plus.bw", sep="")
file.bw.minus.mat <- paste("../",t,"_MB6_all_R1.mat_1bp_minus.bw", sep="")
hl.bw.plus.pat <- paste("../",at,"_MB6_all_R1.pat_1bp_plus.bw", sep="")
hl.bw.minus.pat <- paste("../",at,"_MB6_all_R1.pat_1bp_minus.bw", sep="")
hl.bw.plus.mat <- paste("../",at,"_MB6_all_R1.mat_1bp_plus.bw", sep="")
hl.bw.minus.mat <- paste("../",at,"_MB6_all_R1.mat_1bp_minus.bw", sep="")
heatmap.AT3(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=50000, step=500,  file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT50Kb_step500_navg",navg,".pdf",sep=""), hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat, bl_wd=1, breaks=seq(0,10,1), show.AT.line=TRUE, navg=navg)
heatmap.AT3(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=100000, step=1000, file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT100Kb_step1K_navg",navg,".pdf",sep=""), hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat, bl_wd=1, breaks=seq(0,10,1), show.AT.line=TRUE, navg=navg)
heatmap.AT3(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=200000, step=1000,file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT200Kb_step1K_navg",navg,".pdf",sep=""), hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat, bl_wd=1, breaks=seq(0,10,1), show.AT.line=TRUE, navg=navg)
heatmap.AT3(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=50000, step=500,  file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT50Kb_step500_navg",navg,"_noLine.pdf",sep=""), hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat, bl_wd=1, breaks=seq(0,10,1), show.AT.line=FALSE, navg=navg)
heatmap.AT3(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=100000, step=1000, file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT100Kb_step1K_navg",navg,"_noLine.pdf",sep=""), hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat, bl_wd=1, breaks=seq(0,10,1), show.AT.line=FALSE, navg=navg)
heatmap.AT3(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=200000, step=1000,file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap_AT200Kb_step1K_navg",navg,"_noLine.pdf",sep=""), hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat, bl_wd=1, breaks=seq(0,10,1), show.AT.line=FALSE, navg=navg)
}



# between_organs_AT_bw_heatmap_3("BN","BN")
# between_organs_AT_bw_heatmap_3("LV","LV")
# between_organs_AT_bw_heatmap_3("LV","BN")
# between_organs_AT_bw_heatmap_3("BN","LV")

# between_organs_AT_bw_heatmap_3("BN","BN", navg=10)
# between_organs_AT_bw_heatmap_3("LV","LV", navg=10)
# between_organs_AT_bw_heatmap_3("LV","BN", navg=10)
# between_organs_AT_bw_heatmap_3("BN","LV", navg=10)

# between_organs_AT_bw_heatmap_3("BN","BN", navg=20)
# between_organs_AT_bw_heatmap_3("LV","LV", navg=20)
# between_organs_AT_bw_heatmap_3("LV","BN", navg=20)
# between_organs_AT_bw_heatmap_3("BN","LV", navg=20)

between_organs_AT_bw_heatmap_4 <- function (at, t, navg=1){
    # use Alleleic Termination from one organ and see the proseq reads abundance from another organ
AT_intersect <- read.table(paste("../",at,"_AT_4tunitIntersectNativeHMM_intersectRegion.bed", sep=""), header = F)
file.bw.plus.pat <- paste("../",t,"_all_plus.rpm.bw", sep="")
file.bw.minus.pat <- paste("../",t,"_all_minus.rpm.bw", sep="")
file.bw.plus.mat <- paste("../",t,"_all_plus.rpm.bw", sep="")
file.bw.minus.mat <- paste("../",t,"_all_minus.rpm.bw", sep="")

#heatmap.AT4_allreads(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=50000, step=500,  file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap-allReads_AT50Kb_step500_navg",navg,".pdf",sep=""), bl_wd=1, show.AT.line=TRUE, navg=navg, use.log=TRUE, breaks=seq(0, 5, 1))
heatmap.AT4_allreads(AT_intersect, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat, dist=30000, step=500,file.pdf=paste(at,"-ATtunitAlleleHMMInterectRegion_", t,"-bw-heatmap-allReads_AT30Kb_step500_navg",navg,"_time30_break10.pdf",sep=""), bl_wd=1, show.AT.line=TRUE, navg=navg, times=30, use.log=FALSE, breaks=seq(0, 10, 1))
}



# between_organs_AT_bw_heatmap_4("BN","BN", navg=1)
# between_organs_AT_bw_heatmap_4("LV","LV", navg=1)
# between_organs_AT_bw_heatmap_4("LV","BN", navg=1)
# between_organs_AT_bw_heatmap_4("BN","LV", navg=1)




# file.bw.plus.pat <- paste("../",t,"_all_plus.bw", sep="")
# file.bw.minus.pat <- paste("../",t,"_all_minus.bw", sep="")
# file.bw.plus.mat <- paste("../",t,"_all_plus.bw", sep="")
# file.bw.minus.mat <- paste("../",t,"_all_minus.bw", sep="")

