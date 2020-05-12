library("TmCalculator")
#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/GC_content")


High_Low_GC <- function(df, name){
  pdf(paste(name,"_GC_content.pdf", sep = ""))
  #par(mfrow=c(2,1))
  span=0
  #for (span in 0:1){
    h<- NULL
    for (i in 1:NROW(df)){
      a=s2c(as.character(df$V7[i]))
      t <- NULL
      for (j in 1:(41-span)){
        if (span>0){
          t = c(t, GC(c2s(a[j:(j+span)])))
        }else{
          t = c(t, GC(c2s(c(a[j],"A"))))
        }
      }
      h <- rbind(h,t)
    }
    HighTm = apply(h, 2, FUN = mean) 
    #colMeans(m)
    
    l<- NULL
    for (i in 1:NROW(df)){
      a=s2c(as.character(df$V8[i]))
      t <- NULL
      for (j in 1:(41-span)){
        if (span>0){
          t = c(t, GC(c2s(a[j:(j+span)])))
        }else{
          t = c(t, GC(c2s(c(a[j],"A"))))
        }
      }
      l <- rbind(l,t)
    }
    LowTm = apply(l, 2, FUN = mean) 
    if (span>0){
      plot(-20:(20-span), LowTm, col="blue", xlab="distance to maxTSN", ylab=paste("mean of GC(x + next ",span,"bp)",sep=""), 
           xlim = c(-20,20),
           main=name, type="o", ylim=c(min(HighTm,LowTm),max(HighTm,LowTm)), las=1)
    } else {
      plot(-20:(20-span), LowTm*2, col="blue", xlab="distance to maxTSN", ylab="mean of GC content", 
           xlim = c(-20,20),
           main=name, type="o", ylim=c(0,max(HighTm,LowTm)*2), las=1)
      
    }
    points(-20:(20-span), HighTm*2,col="red", type="o")
    legend("bottomleft", legend=c("High", "Low"),
           col=c("red", "blue"), bty = "n", lty=1, pch=1)
  #}
  dev.off()
}

for(t in c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")){
  df=read.table(file =paste(t,"_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1_+-20_High_LowAlleleSeq.bed", sep=""))
  #V7 Allele with high expression
  #V8 Allele with low expression
  for (i in 1:NROW(df)){
    df$HighAlleleGC[i] = GC(as.character(df$V7[i]))
    df$LowAlleleGC[i] = GC(as.character(df$V8[i]))
  }
  pdf(paste(t,"_GC_content_hist.pdf", sep = ""))
  hist(df$HighAlleleGC-df$LowAlleleGC, main=t, col="blue")
 dev.off()
    High_Low_GC(df, name=t)
}