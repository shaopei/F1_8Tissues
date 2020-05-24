#R --vanilla --slave --args BN_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs_+-200_High_LowAlleleSeq.bed BN 10 < getGC_content_HighLowAllele.R &
args=(commandArgs(TRUE))
df_fp=args[1]
t=args[2]
step=as.integer(args[3])

library("TmCalculator")
#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/GC_content")


High_Low_GC <- function(df, name, step=5){
  pdf(paste(name,"_GC_content_step=",step,".pdf", sep = ""))
  #par(mfrow=c(4,1))
  #for (step in c(1,5,10,20)){
  h<- NULL
  for (i in 1:NROW(df)){
    a=s2c(as.character(df$HighAlleleSeq[i]))
    len_a = length(a)
    temp <- NULL
    if (step>1){
      for (j in seq(1,(len_a-step), step)){
        temp = c(temp, GC(c2s(a[j:(j+step-1)])))
      }
      h <- rbind(h,temp)
    }else{
      for (j in seq(1,(len_a), step)){
        temp = c(temp, GC(c2s(c(a[j],"A"))))
      }
      h <- rbind(h,temp)
    }
  }
  HighTm = apply(h, 2, FUN = mean, na.rm = TRUE) 
    #colMeans(m)
    
    l<- NULL
    for (i in 1:NROW(df)){
      a=s2c(as.character(df$LowAlleleSeq[i]))
      temp <- NULL
      if (step>1){
        for (j in seq(1,(len_a-step), step)){
          temp = c(temp, GC(c2s(a[j:(j+step-1)])))
        }
        l <- rbind(l,temp)
      }else{
        for (j in seq(1,(len_a), step)){
          temp = c(temp, GC(c2s(c(a[j],"A"))))
        }
        l <- rbind(l,temp)
      }
    }

    LowTm = apply(l, 2, FUN = mean, na.rm = TRUE) 
    show.window = floor(len_a/2)
    if (step >1){
      x <- seq(-1*(show.window)+floor(step/2), show.window, step)
    }else{
      x <- seq(-1*(show.window), show.window, step)
    }
    
    if (step>1){
      plot(x[1:length(LowTm)], LowTm, col="blue", xlab="distance to as.maxTSN", ylab=paste("mean of GC(around x, bin=",step,"bp)",sep=""),
           #xlim=c(-50,50),
           main=name, type="o", ylim=c(min(HighTm,LowTm),max(HighTm,LowTm)), las=1)
      points(x[1:length(LowTm)], HighTm,col="red", type="o")
      } else {
      plot(x[1:length(LowTm)], LowTm*2, col="blue", xlab="distance to maxTSN", ylab="mean of GC content", 
           #xlim=c(-50,50),
           main=name, type="o", ylim=c(min(HighTm*2,LowTm*2),max(HighTm*2,LowTm*2)), las=1)
      points(x[1:length(LowTm)], HighTm*2,col="red", type="o")
    }
   
    legend("bottomleft", legend=c("High", "Low"),
           col=c("red", "blue"), bty = "n", lty=1, pch=1)
  #}
  dev.off()
  return (list(x, HighTm, LowTm))
}

if(0){
for(t in c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")){
  df=read.table(file =paste(t,"_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr0.1_+-20_High_LowAlleleSeq.bed", sep=""))
  colnames(df)[7:8] = c("HighAlleleSeq", "LowAlleleSeq")
  #V7 Allele with high expression
  #V8 Allele with low expression
  for (i in 1:NROW(df)){
    df$HighAlleleGC[i] = GC(as.character(df$HighAlleleSeq[i]))
    df$LowAlleleGC[i] = GC(as.character(df$LowAlleleSeq[i]))
  }
  pdf(paste(t,"_GC_content_hist.pdf", sep = ""))
  hist(df$HighAlleleGC-df$LowAlleleGC, main=t, col="blue")
 dev.off()
    High_Low_GC(df, name=t)
}
}

#for(t in c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")){
  name=t
  #df=read.table(file =paste(t,"_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs_+-200_High_LowAlleleSeq.bed", sep=""))
  df=read.table(file =df_fp)

  colnames(df)[16:17] = c("HighAlleleSeq", "LowAlleleSeq")
  colnames(df)[7:10]=c("masked_p_value","masked_fdr", "unmasked_p_value","unmasked_fdr")
  colnames(df)[11] = "winP"
  # for (i in 1:NROW(df)){
  #   df$HighAlleleGC[i] = GC(as.character(df$HighAlleleSeq[i]))
  #   df$LowAlleleGC[i] = GC(as.character(df$LowAlleleSeq[i]))
  # }
  # pdf(paste(t,"_GC_content_hist.pdf", sep = ""))
  # x=df$HighAlleleGC-df$LowAlleleGC
  # hist(x, main=t, col="blue", breaks = seq(-3.005,3,0.01), xlab="df$HighAlleleGC-df$LowAlleleGC")
  # dev.off()
  
  ###
  g9 = df[df$unmasked_fdr> 0.9,]
  p_value_cut = max(df$unmasked_p_value[df$unmasked_fdr<=0.1])
  
  # asTSS driven by multiple base
  m = df[df$masked_p_value <= p_value_cut & df$unmasked_fdr<=0.1,]
  # asTSS driven by single base
  s = df[df$masked_p_value > p_value_cut & df$unmasked_fdr<=0.1,]
  cat(name, "\n")
  cat("dim(g9); max(g9$unmasked_p_value); min(g9$unmasked_p_value)", "\n")
  cat (dim(g9), max(g9$unmasked_p_value), min(g9$unmasked_p_value) )
  cat ("dim(m), max(m$masked_p_value), min(m$masked_p_value)", "\n")
  cat (dim(m), max(m$masked_p_value), min(m$masked_p_value), "\n")
  cat("  dim(s); max(s$masked_p_value); min(s$masked_p_value)", "\n")
  cat(dim(s), max(s$masked_p_value), min(s$masked_p_value), "\n")
  
 # for (step in c(1,5, 9,13,19)){
    g9_gc=High_Low_GC(g9, name=paste(t,"_fdr>0.9_step=",step,sep=""), step=step)
    m_gc=High_Low_GC(m, name=paste(t,"_m_step=",step,sep=""), step=step)
    s_gc=High_Low_GC(s, name=paste(t,"_s_step=",step,sep=""), step=step)
    
    pdf(paste(t,"_deltaGC_step=",step,".pdf" ,sep=""), width = 10, height = 10)
    plot(s_gc[[1]][1:length(s_gc[[2]])], s_gc[[2]]-s_gc[[3]], type="o", col="red", 
         ylim=c(min(s_gc[[2]]-s_gc[[3]], m_gc[[2]]-m_gc[[3]]), max(s_gc[[2]]-s_gc[[3]], m_gc[[2]]-m_gc[[3]])),
         main=paste(t," Mean GC in High Allele - Mean GC Low Allele, bin size=",step,sep=""),
         xlab="distance to maxTSN", ylab="delta GC", pch=19)
    abline(h=0)
    abline(v=0)
    points(g9_gc[[1]][1:length(g9_gc[[2]])], g9_gc[[2]]-g9_gc[[3]], type="o", col="gray", pch=19)
        points(m_gc[[1]][1:length(m_gc[[2]])], m_gc[[2]]-m_gc[[3]], type="o", col="blue", pch=19)

    points(s_gc[[1]][1:length(s_gc[[2]])], s_gc[[2]]-s_gc[[3]], type="o", col="red", pch=19)
    
    
    legend("bottomleft", legend=c("fdr > 0.9","fdr <= 0.1, Single", "fdr <= 0.1, Multiple"),
           col=c("gray","red", "blue"), bty = "n", lty=1, pch=19)
    dev.off()
  #}
  
#}


#df$deltaGC = df$HighAlleleGC-df$LowAlleleGC
