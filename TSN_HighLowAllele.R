library("TmCalculator")
library(bigWig)

file_dir="~/Box Sync/BN_IGV/"
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/TSN_ShootingGallery/")

getTSN_2N <- function(seq, d){
  seq=s2c(as.character(seq))
  return (c2s(seq[d:(d+1)]))
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

organ="BN"; d=50; name=organ
Dist <- NULL
Delta_Signal <- NULL
for (strand in c("+", "-")){
  #df=read.table(paste(organ, "_allReads_TSN5+_SNP_binomtest_interestingHets_+-",d,"_High_LowAlleleSeq.bed", sep=""))
  df=read.table(paste(organ, "_allReads_TSN5+_SNP_binomtest_+-",d,"_High_LowAlleleSeq.bed", sep=""))
  colnames(df)[12:13] = c("HighAlleleSeq", "LowAlleleSeq")
  colnames(df)[9]="Bino_p_value"
  colnames(df)[11] = "winP"
  
  
  
  
  df$CAatHighTSN =  sapply(df$HighAlleleSeq, getTSN_2N, d=d)
  df$CAatLowTSN =  sapply(df$LowAlleleSeq, getTSN_2N, d=d)
  
  # keep High allele with CA at TSN
  df = df[df$CAatHighTSN=="CA",]
  df = df[df$V10 == strand,]
  
  # keep rows with total allelic reads > 5
  df = df[df$V6+df$V7 >= 5,]
  
  # make a matrix of dinucleotide. Sliding windows with overlaps
  h_m <- NULL
  for (i in 1:length(df$HighAlleleSeq)){
    a = s2c(as.character(df$HighAlleleSeq[i]))
    r <- NULL
    for (j in seq(1, 2*d, 1)){
      r = c(r, c2s(a[j:(j+1)]))
    }
    h_m  = rbind(h_m ,r)
  }
  
  l_m <- NULL
  for (i in 1:length(df$LowAlleleSeq)){
    a = s2c(as.character(df$LowAlleleSeq[i]))
    r <- NULL
    for (j in seq(1, 2*d, 1)){
      r = c(r, c2s(a[j:(j+1)]))
    }
    l_m  = rbind(l_m ,r)
  }
  
  dim(h_m)
  dim(l_m)
  
  diNu = c("CA", "CG", "TA", "TG")
  

  
  # the transcription level neat the TSN
  TSS = df[,c(1,2,3,4,5,10)]
  if (strand=="+"){
    TSS[,2] <- df$V2 - d + 1
    TSS[,3] <- df$V3 + d
  }else{
    TSS[,2] <- df$V2 - d 
    TSS[,3] <- df$V3 + d - 1
  }
  
  step=1; times=1; use.log=FALSE
  
  file.bw.plus.pat=paste(file_dir,"BN_map2ref_1bpbed_map5_CAST_plus.bw", sep="")
  file.bw.minus.pat=paste(file_dir,"BN_map2ref_1bpbed_map5_CAST_minus.bw", sep="")
  file.bw.plus.mat=paste(file_dir,"BN_map2ref_1bpbed_map5_B6_plus.bw", sep="")
  file.bw.minus.mat=paste(file_dir,"BN_map2ref_1bpbed_map5_B6_minus.bw", sep="")
  readCount.pat <- read_read_mat_S (file.bw.plus.pat, file.bw.minus.pat, TSS[,c(1:6)], step, times=times, use.log=use.log)
  readCount.mat <- read_read_mat_S (file.bw.plus.mat, file.bw.minus.mat, TSS[,c(1:6)],  step, times=times, use.log=use.log) 
  readCount.combined <- readCount.pat + readCount.mat
  dim(readCount.pat)
  
  # generate x axis of the plot, distance to the TSN
  x_m  <- NULL
  for (i in 1:dim(readCount.pat)[1]){
    x_m <- rbind( x_m, seq(-d,d-1,1))
  }
  dim(x_m)
  t=1
  # delta_reads =  log2((readCount.HIGH+t)/(readCount.LOW+t))
  delta_reads = log2((readCount.pat+t)/(readCount.mat+t))
  delta_reads[which(df$winP=="M"),] = log2((readCount.mat+t )/( readCount.pat+t))[which(df$winP=="M"),]
  #delta_reads = ((readCount.pat+t)-(readCount.mat+t))
  #delta_reads[which(df$winP=="M"),] = ((readCount.mat+t )-( readCount.pat+t))[which(df$winP=="M"),]
  dim(delta_reads)
  #View(delta_reads)
  #dim(df)

  
  # matrix indicates if the dinucleotide is "CA", "CG", "TA", or "TG"
  m = matrix(h_m %in% diNu & readCount.combined >=5, ncol = 100)
  dim(m)
  
  
  if (strand=="+"){
    plot(x_m[m], delta_reads[m], col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19, 
         #ylim=c(-100,100),
         ylab = "log2(High Allele + 1 / Low Allele + 1)",
         xlab = "Distance to TSN with CA at High Allele")
    Dist <- c(Dist,x_m[m] )
    Delta_Signal <- c(Delta_Signal, delta_reads[m])
    
  }else{
    points(x_m[m], delta_reads[m], col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19)
    Dist <- c(Dist, x_m[m])
    Delta_Signal <- c(Delta_Signal, delta_reads[m])
  }
}

out = data.frame(dist=Dist, y= Delta_Signal)
out = out[out$dist!=-1,]
signal.lo <-  loess(y ~ dist, out)
x=seq(-d,d-1,1)
p=predict(signal.lo, data.frame(dist = x), se = TRUE)
plot(out, col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.15), pch=19,
     ylab = "log2(High Allele + 1 / Low Allele + 1)",
     xlab = "Distance to TSN with CA at High Allele")#, ylim=c(5,-5))
lines(x, p$fit, col="red" )
plot(x, p$fit, col="red",  ylab = "Predicted log2(High Allele + 1 / Low Allele + 1)",
     xlab = "Distance to TSN with CA at High Allele")
abline(h=0)

  
  plot(x_m, abs(readCount.pat-readCount.mat), col=rgb(red=0.2, green=0.2, blue=0.2, alpha=0.2), pch=19)#, ylim=c(0,150))
  