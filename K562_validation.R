setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/K562_validation/")

read_read_mat_S <-function (file.plus.bw, file.minus.bw , bed6, step=1, navg = 1, times=1, use.log=FALSE)
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
  

metaplot <-function (bed_file){
  bed6=read.table(bed_file)
  bed6=bed6[bed6$V1 !="chrM",]
  AT <- bed6
  d=30
  for (i in 1:NROW(bed6)){
    if(bed6[i,6]=="-") {
      bed6[i,3] <- AT[i,3] + d
      bed6[i,2] <- AT[i,3] - d}
    else
    {bed6[i,2] <- AT[i,2] - d
    bed6[i,3] <- AT[i,2] + d
    }
  }
  #bed6$V3 - bed6$V2
  file.bw.plus = "groseq_tss_wTAP_plus.bigWig"
  file.bw.minus = "groseq_tss_wTAP_minus.bigWig"
  step=1
  times=1
  use.log=FALSE
  hmat.high <- read_read_mat_S (file.bw.plus, file.bw.minus, bed6[,c(1:6)], step, times=times, use.log=use.log)
  
  if(NROW(bed6)<100){
    meta.hmat.high <- colMedians(hmat.high)
  }else{
    meta.hmat.high.temp <- subsampled.quantiles.metaprofile(hmat.high)
    meta.hmat.high <- meta.hmat.high.temp$middle
  }
  
  show.window=d
  return(meta.hmat.high)
}

pos <- metaplot("promoterNenhancers.maxTSS.K562.bed")
tss <- metaplot("ChROseq_merged_0h_TSS_maxTSNs.bed")
neg <- metaplot("gencode.V19.annotation_5END.bed")


pdf("K562_metaplot.pdf", width=10, height = 8, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)


plot(seq(-29.5,30,1), pos, type="l", ylab="GRO-cap signal", xlab="distance to maxTSS", las=1)
lines(seq(-29.5,30,1), tss, col="blue")
lines(seq(-29.5,30,1), neg, col="red")
legend("topright", 
       legend = c("Co-PRO", "TSS", "gencode"),
       #pch=c(15,15),
       cex=2, 
       lty=1,
       #bty="n",
       lwd=1.5, 

       #angle=45,
       col=c("black","blue", "red")
       , bty = "n"
)

genecode_meta.hmat.high = meta.hmat.high 

