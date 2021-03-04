require(bigWig)
require(pheatmap)
require(RColorBrewer)
library(gtable)
library(grid)
library(sqldf)
library("vioplot")

# functions

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


heatmap.AT3 <-function(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat ,
                       dist, step, up_dist =20000, file.pdf="heatmap.pdf", 
                       hl.bw.plus.pat=NULL,hl.bw.minus.pat=NULL, hl.bw.plus.mat=NULL ,hl.bw.minus.mat=NULL, 
                       bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 10, use.log=FALSE, times=1){
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
    # based on file.bw.plus.pat, file.bw.minus.pat, file.bw.plus.mat, file.bw.minus.mat
    hmat.pat.AT.rowSums <- NULL
    hmat.mat.AT.rowSums <- NULL
    for (i in 1:NROW(hmat.pat)){
      hmat.pat.AT.rowSums[i] <- sum(hmat.pat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
      hmat.mat.AT.rowSums[i] <- sum(hmat.mat[i,][AT$start.steps[i]:min(AT$end.steps[i],bin_number)]) 
    }
  } else {
    # dertemine High/Low allele based on the reads count within the AT regions (Before adjust by dist)
    # based on hl.bw.plus.pat, hl.bw.minus.pat, hl.bw.plus.mat, hl.bw.minus.mat
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

heatmap.AT4_allreads <-function(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat ,
                                dist, step, up_dist =20000, file.pdf="heatmap.pdf", 
                                hl.bw.plus.pat=NULL,hl.bw.minus.pat=NULL, hl.bw.plus.mat=NULL ,hl.bw.minus.mat=NULL,
                                bl_wd=1, breaks=NULL, cols=NULL, show.AT.line=TRUE, navg = 10, use.log=FALSE, times=1){
  # file.bw.plus.pat = file.bw.plus.mat AND , file.bw.minus.pat= file.bw.minus.mat
  # used to get heatmap from all reads
  
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

read_count <- function(file.plus.bw, file.minus.bw , bed6) {
  bw.plus  <- load.bigWig( file.plus.bw )
  bw.minus <- load.bigWig( file.minus.bw )
  CountMatrix <- bed6.region.bpQuery.bigWig(bw.plus, bw.minus, bed6[,c(1:6)] ,abs.value=TRUE, op = "sum")
  unload.bigWig(bw.plus);
  unload.bigWig(bw.minus);
  
  return(CountMatrix);    
}



Tversion <- "4tunitIntersectNativeHMM"
getTunitLength <- function(t, path="./", body=paste("_AT",Tversion, "tunits.bed", sep = "_")){
  df <- read.table(paste(path,t,body, sep=""), header = F)
  return(log10(df$V3-df$V2))
}

getATLength <- function(t, path="./", body=paste("_AT",Tversion, "intersectRegion.bed", sep = "_")){
  df <- read.table(paste(path,t,body, sep=""), header = F)
  return(log10(df$V3-df$V2))
} 


plot_colorByDensity = function(x1,x2,
                               ylim=c(min(x2),max(x2)),
                               xlim=c(min(x1),max(x1)),
                               xlab="",ylab="",main="", bty="n") {
  
  df <- data.frame(x1,x2)
  x <- densCols(x1,x2, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  cols <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
  df$col <- cols[df$dens]
  plot(x2~x1, data=df[order(df$dens),],
       ylim=ylim,xlim=xlim,pch=20,col=col,
       cex=2,xlab=xlab,ylab=ylab,
       main=main, bty=bty, las=1)
}

setwd(paste("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/gencode.vM25/tunitIntersectNativeHMM",Tversion,sep = "/"))
tissue_list_full <-c("BN","HT", "SK", "LV","GI", "ST" , "KD", "SP")

### the relationship between Tunit expression level and AT window length
pdf("scatter_plot_logATRPKM_logTunitRPKM.pdf",width=16, height = 8, useDingbats=FALSE)
par(mfrow=c(2,4))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2, cex.main=2.5)
for (at in tissue_list_full){  
  #for (at in c("BN",  "SP", "LV")){
  tunit <-  read.table(paste("./",at,"_AT_",Tversion, "_tunits.bed", sep=""), header = F)
  AT <- read.table(paste("./",at,"_AT_",Tversion,"_intersectRegion.bed", sep=""), header = F)
  file.plus.bw <- paste("~/Box Sync/IGV/bigWig_notBN/",at,"_all_plus.rpm.bw", sep="")
  file.minus.bw <- paste("~/Box Sync/IGV/bigWig_notBN/",at,"_all_minus.rpm.bw", sep="")
  tunit_rpkm <- read_count(file.plus.bw, file.minus.bw, tunit) *1000/(tunit$V3 - tunit$V2) 
  AT_rpkm <- read_count(file.plus.bw, file.minus.bw, AT) *1000/(AT$V3 - AT$V2) 
  
#  plot(log10(tunit_rpkm), log10(AT_rpkm), main=paste(at,Tversion,sep=" "), xlim = c(-2,2), ylim = c(-2,2))

  plot_colorByDensity(log10(tunit_rpkm), log10(AT_rpkm), main=at,
                      xlim = c(-2,2), ylim = c(-2,2),
                      xlab="log10(tunit RPKM)", ylab=" log10(AT RPKM)")
}
dev.off()

# heat map

navg=1

tissues <- c("BN", "LV","SP")

Tunit_BW_pairs <- function (at, t){
  tunit <-  read.table(paste("./",at,"_AT_",Tversion, "_tunits.bed", sep=""), header = F)
  file.plus.bw <- paste("~/Box Sync/IGV/bigWig_notBN/",t,"_all_plus.rpm.bw", sep="")
  file.minus.bw <- paste("~/Box Sync/IGV/bigWig_notBN/",t,"_all_minus.rpm.bw", sep="")
  log10(read_count(file.plus.bw, file.minus.bw, tunit)*1000/(tunit$V3 - tunit$V2)+0.01) # log10(rpkm+0.01)
}

# Fig 5B
at="LV"
for (t in tissues){
  AT <- read.table(paste("./",at,"_AT_",Tversion,"_intersectRegion.bed", sep=""), header = F)
  t.value <-Tunit_BW_pairs (at, t)
  
  # use MB6 and PB6 (n=7) BN_map2ref_1bpbed_map5_B6_minus.bw
  file.bw.plus.pat  <- paste("~/Box Sync/BN_IGV/",t,"_map2ref_1bpbed_map5_CAST_plus.bw", sep="")
  file.bw.minus.pat <- paste("~/Box Sync/BN_IGV/",t,"_map2ref_1bpbed_map5_CAST_minus.bw", sep="")
  file.bw.plus.mat  <- paste("~/Box Sync/BN_IGV/",t,"_map2ref_1bpbed_map5_B6_plus.bw", sep="")
  file.bw.minus.mat <- paste("~/Box Sync/BN_IGV/",t,"_map2ref_1bpbed_map5_B6_minus.bw", sep="")
  hl.bw.plus.pat    <- paste("~/Box Sync/BN_IGV/",at,"_map2ref_1bpbed_map5_CAST_plus.bw", sep="")
  hl.bw.minus.pat   <- paste("~/Box Sync/BN_IGV/",at,"_map2ref_1bpbed_map5_CAST_minus.bw", sep="")
  hl.bw.plus.mat    <- paste("~/Box Sync/BN_IGV/",at,"_map2ref_1bpbed_map5_B6_plus.bw", sep="")
  hl.bw.minus.mat   <- paste("~/Box Sync/BN_IGV/",at,"_map2ref_1bpbed_map5_B6_minus.bw", sep="")
  file.plus.bw      <- paste("~/Box Sync/BN_IGV/",at,"_map2ref_1bpbed_map5_plus.bw", sep="")
  file.minus.bw     <- paste("~/Box Sync/BN_IGV/",at,"_map2ref_1bpbed_map5_minus.bw", sep="")
  
  
  # expression of all reads (not only allelic specific ones)
  heatmap.AT4_allreads(AT, file.plus.bw,file.minus.bw, file.plus.bw ,file.minus.bw, 
                       dist=20000, step=500,
                       file.pdf=paste(at,".bw_",1,".pdf",sep=""),
                       bl_wd=1, show.AT.line=TRUE, navg=navg, times=1, use.log=FALSE, breaks=seq(0, 20, 0.01))
  # expression of allelic specific reads
  # only use tunits with  log10(rpkm+0.01) > -1
  # high low allele based on at (the ChRO-seq reads)
  heatmap.AT3(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat,
              hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat,
              dist=20000, step=500,
              file.pdf=paste(at,".AT_",t,".bw_",2,".pdf",sep=""),
              bl_wd=1, show.AT.line=TRUE, navg=navg, use.log=FALSE,
              times=10, breaks = seq(0,20,0.01))
  heatmap.AT3(AT, file.bw.plus.pat,file.bw.minus.pat, file.bw.plus.mat ,file.bw.minus.mat,
              hl.bw.plus.pat=hl.bw.plus.pat ,hl.bw.minus.pat=hl.bw.minus.pat, hl.bw.plus.mat = hl.bw.plus.mat ,hl.bw.minus.mat=hl.bw.minus.mat,
              dist=20000, step=500,
              file.pdf=paste(at,".AT_",t,".bw_",2,"_noline.pdf",sep=""),
              bl_wd=1, show.AT.line=FALSE, navg=navg, use.log=FALSE,
              times=10, breaks = seq(0,20,0.01))
}

# Fig 5C
# get the length distribution of AT windows from all organs
f1_p = "T8_AT_4tunitIntersectNativeHMM_intersectRegion.bed"
pdf("T8_ATwindow_length.pdf", width=5, height = 5, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
f1=read.table(f1_p,header=F)
h<- hist(log10(f1$V3-f1$V2)
         , breaks = seq(0,7,1)
         ,plot =FALSE
         ,right = FALSE)

h$counts=h$counts/sum(h$counts)
plot(h,
     col="red", 
     #,density=25     
     ylab="Proportion",
     xlim=c(0,7),
     xlab="AT window length (bp)",
     las=2,
     xaxt='n',
     main=""
     )
axis(1, at=seq(0,7,1), labels=c(0,"10","102","103","104","105","106","107"), las=1)
# axis(1, at=seq(0,7,1), labels=c(0,10,100,"1,000","10,000","100,000","1,000,000","10,000,000"), las=2)

dev.off()

# Sup Fig 5B
tussue_list <-c("BN","HT", "SK", "LV","GI", "ST" , "KD", "SP")
pdf("Organs_ATwindow_length_violin.pdf", width=10, height = 8, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)

  getATLength <- function(t, path="./", body=paste("_AT",Tversion, "intersectRegion.bed", sep = "_")){
    df <- read.table(paste(path,t,body, sep=""), header = F)
    return(log10(df$V3-df$V2))
  } 
  
  vioplot(getATLength(tussue_list[1]),getATLength(tussue_list[2]),getATLength(tussue_list[3]),getATLength(tussue_list[4]),getATLength(tussue_list[5]),getATLength(tussue_list[6]),getATLength(tussue_list[7]),getATLength(tussue_list[8]),
          names = tussue_list ,
          ylab= "log10(AT window length)",
#          ylim = c(0,6),
          main="" ,
          col="red",
          las=1)

  dev.off()
  