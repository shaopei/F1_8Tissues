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


read_read_mat_S <-function (file.plus.bw, file.minus.bw , bed6, step=1, navg = 20, times=1, use.log=FALSE)
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

pause_window <- read.table("/Volumes/SPC_SD/KD_IGV/KD_dREG_5mat5pat_uniq_pValue_fdr0.1.bed.bed", header = F)
AT <- pause_window
AT <- AT[AT$V1 != 'chrX',]
t="KD"
#KD_PB6_F5_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_plus.bw
file.bw.plus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_plus.bw", sep="")
file.bw.minus.pat <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5_dedup_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_minus.bw", sep="")
file.bw.plus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_plus.bw", sep="")
file.bw.minus.mat <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5_dedup_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.1bp.sorted_minus.bw", sep="")
file.plus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5_dedup_QC_end_plus.bw", sep="")
file.minus.bw <- paste("/Volumes/SPC_SD/KD_IGV/",t,"_PB6_F5_dedup_QC_end_minus.bw", sep="")


breaks =seq(0, 10, 0.001)

heatmap.Pause <-function(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat ,dist=1000, step=2, up_dist =100, file.pdf="heatmap.pdf", hl.bw.plus.pat=NULL,hl.bw.minus.pat=NULL, hl.bw.plus.mat=NULL ,hl.bw.minus.mat=NULL, bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 10, use.log=FALSE, times=1){
  # AT here is the window of pause
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
  hmat.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, bed6[,c(1:6)], step, times=times, use.log=use.log)
  hmat.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, bed6[,c(1:6)],  step, times=times, use.log=use.log)     
  
  
  
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





heatmap.Pause_allreads <-function(AT, file.plus.bw, file.minus.bw ,dist=1000, step=2, up_dist =100, file.pdf="heatmap.pdf", hl.bw.plus.pat=NULL,hl.bw.minus.pat=NULL, hl.bw.plus.mat=NULL ,hl.bw.minus.mat=NULL, bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 10, use.log=FALSE, times=1){
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

