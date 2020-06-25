
library("TmCalculator")
library(seqLogo)
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/GC_content")



#define function that divides the frequency by the row sum i.e. proportions
proportion <- function(x){
  rs <- sum(x);
  return(x / rs);
}

seq_upperCase <- function(seq){
  seq[seq=="a"]<- "A"
  seq[seq=="t"]<- "T"
  seq[seq=="c"]<- "C"
  seq[seq=="g"]<- "G"
  
  a <- NULL
  t <- NULL
  c <- NULL
  g <- NULL
  for (i in 1:NCOL(seq)){
    a <- c(a, sum(seq[,i]=="A"))
    t <- c(t, sum(seq[,i]=="T"))
    c <- c(c, sum(seq[,i]=="C"))
    g <- c(g, sum(seq[,i]=="G"))
  }
  return (data.frame(a,c,g,t))
}

SeqLogo <- function(seq, output, range) {
  #seq=m$HighAlleleSeq
  seq<- data.frame(do.call(rbind, strsplit(as.character(seq), "")))
  df <- seq_upperCase(seq)
  #create position weight matrix
  pwm <- apply(df, 1, proportion)
  p = makePWM((pwm[,range]))
  #p <- makePWM(pwm)
  # slotNames(p)
  # p@consensus
  # p@ic
  # p@width
  # p@alphabet
  pdf(output)
  seqLogo(p)
  dev.off()
  return (pwm)
}

High_Low_GC_excludeCA <- function(df, name, step=5){
  pdf(paste(name,"_GC_content_excludeCA_step=",step,".pdf", sep = ""))
  #par(mfrow=c(4,1))
  #for (step in c(1,5,10,20)){
  h<- NULL
  for (i in 1:NROW(df)){
    a=s2c(as.character(df$HighAlleleSeq[i]))
    len_a = length(a)
    d = (len_a-1)/2
    #CA at a[d:(d+1)]
    a = c(a[1:(d-1)] ,"N","N",a[(d+2):len_a])
    #len_a = length(a)
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
    len_a = length(a)
    d = (len_a-1)/2
    #CA at a[d:(d+1)]
    a = c(a[1:(d-1)] ,"N","N",a[(d+2):len_a])
    #len_a = length(a)
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
    plot(x[1:length(LowTm)], LowTm, col="blue", xlab="distance to as.maxTSN", 
         ylab=paste("mean of GC(around x, bin=",step,"bp)",sep=""),
         #xlim=c(-50,50),
         main=paste(name, " n=", dim(df)[1], sep=""), 
         type="o", 
         ylim=c(min(HighTm,LowTm),max(HighTm,LowTm)), 
         las=1)
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
  return (list(x, HighTm, LowTm, show.window))
}

df_process <- function(df, name, step){
  
  colnames(df)[16:17] = c("HighAlleleSeq", "LowAlleleSeq")
  colnames(df)[7:10]=c("masked_p_value","masked_fdr", "unmasked_p_value","unmasked_fdr")
  colnames(df)[11] = "winP"
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
  
  
  ### seglogo
  
  cat("start SeqLogo\n")
  #d=s_gc[[4]]
  w=10
  range=(d+1-w): (d+1+w)
  organ=name
  m_h = SeqLogo(m$HighAlleleSeq, paste(organ,"_m_HighAllele.pdf", sep=""),range)
  m_l = SeqLogo(m$LowAlleleSeq, paste(organ,"_m_LowAllele.pdf", sep=""),range)
  s_h = SeqLogo(s$HighAlleleSeq, paste(organ,"_s_HighAllele.pdf", sep=""),range)
  s_l = SeqLogo(s$LowAlleleSeq, paste(organ,"_s_LowAllele.pdf", sep=""),range)
  g9_h = SeqLogo(g9$HighAlleleSeq, paste(organ,"_g9_HighAllele.pdf", sep=""),range)
  g9_l = SeqLogo(g9$LowAlleleSeq, paste(organ,"_g9_LowAllele.pdf", sep=""),range)   
  
  acgt_col=c("dark green", "blue", "orange" , "red")
  acgt=c("A","C","G","T")
  
  pdf(paste(organ,"_deltaATCG_single.pdf" ,sep=""), width =7, height = 7)
  #par(mfcol=c(2,1))
  par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  par(mgp=c(3,1,0))
  par(cex.lab=2.2, cex.axis=2.2)
  #par(cex=1.5)
  pch_u=15
  plot(-w:w,s_h[1,range] - s_l[1,range], col = acgt_col[1], type="o",
       ylim=c(-0.1,0.1), pch=pch_u,
       ylab="HighAllele - LowAllele",
       xlab="Distance to maxTSN",
       main=paste(organ,"fdr<=0.1 single", sep=" "),
       las=1, frame=FALSE
  )
  abline(h=0, col="gray")
  for (i in 4:1){
    points(-w:w,s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o", pch=pch_u)
  }
  legend("topleft", legend=acgt,
         col=acgt_col, bty = "n", lty=1, pch=pch_u)
  
  if(0){
    plot(-w:w,m_h[1,range] - m_l[1,range], col = acgt_col[1], type="o",
         ylim=c(-0.1,0.1), pch=1,
         ylab="HighAllele - LowAllele",
         xlab="Distance to maxTSN",
         main=paste(organ, "fdr<=0.1 multiple", sep=" ")
    )
    abline(h=0, col="gray")
    for (i in 1:4){
      points(-w:w,m_h[i,range] - m_l[i,range], col = acgt_col[i], type="o", pch=1)
    }
    legend("bottomleft", legend=acgt,
           col=acgt_col, bty = "n", lty=1, pch=1)
  }
  
  dev.off()
  
  pdf(paste(organ,"_deltaATCG_s_vs_m-4.pdf" ,sep=""), width = 10, height = 10)
  par(mfcol=c(4,1))
  for (i in 1:4){
    plot(-w:w,s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o",
         ylim=c(min(s_h[i,range] - s_l[i,range], m_h[i,range] - m_l[i,range]),max(s_h[i,range] - s_l[i,range], m_h[i,range] - m_l[i,range])), 
         pch=19,
         ylab="HighAllele - LowAllele",
         xlab="Distance to maxTSN",
         main=paste(organ, acgt[i], sep=" ")
    )
    abline(h=0, col="gray")
    points(-w:w,m_h[i,range] - m_l[i,range], col = "black", type="o", pch=1)
    legend("bottomleft", legend=c("single", "multiple"),
           col=c(acgt_col[i], "black"), bty = "n", lty=1, pch=c(19,1))
  }
  dev.off()
  
  
  
  
  g9_gc_e=High_Low_GC_excludeCA(g9, name=paste(organ,"_fdr>0.9_d=",d,sep=""), step=step)
  m_gc_e=High_Low_GC_excludeCA(m, name=paste(organ,"_m_d=",d,sep=""), step=step)
  s_gc_e=High_Low_GC_excludeCA(s, name=paste(organ,"_s_d=",d,sep=""), step=step)
  return(list(g9_gc_e, m_gc_e, s_gc_e, dim(g9)[1], dim(m)[1], dim(s)[1] ))
}

d=35; step=5

organ="BN"
BN_df=read.table(paste(organ, "_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs_+-",d,"_High_LowAlleleSeq.bed", sep=""))
#name=organ

BN_output = df_process(BN_df, organ, step)

organ="LV"; 
LV_df=read.table(paste(organ, "_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs_+-",d,"_High_LowAlleleSeq.bed", sep=""))
LV_output = df_process(LV_df, organ, step)



organ="BN"
output = BN_output
g9_gc_e = output[[1]]
m_gc_e = output[[2]]
s_gc_e = output[[3]]
n_g9 = output[[4]] 
n_m = output[[5]] 
n_s = output[[6]]

### Single base plot
# width=9.23, height = 6.52)
pdf(paste(organ,"_deltaGC_exclude_CA_step=",step,"_d=",d,".pdf" ,sep=""), width=9.23, height = 6.52/2)
#par(mfcol=c(2,1))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
#par(cex=1.5)
pch_u=15
#Brain Single base plot
plot(s_gc_e[[1]][1:length(s_gc_e[[2]])], s_gc_e[[2]]-s_gc_e[[3]], type="o", col="red", 
     #     ylim=c(min(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]]), max(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]])),
     main=paste(organ, " High Allele - Low Allele, d=",g9_gc_e[[4]] ," bin size=",step,"_exclude maxTSN",sep=""),
     xlab="distance to maxTSN", ylab="delta GC", pch=pch_u, xlim=c(-d,d), las=1, frame=FALSE)
#abline(h=0)
#abline(v=0)



points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], g9_gc_e[[2]]-g9_gc_e[[3]], type="o", col="black", pch=pch_u)
#points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=19)

points(s_gc_e[[1]][1:length(s_gc_e[[2]])], s_gc_e[[2]]-s_gc_e[[3]], type="o", col="red", pch=pch_u)

legend("bottomleft", legend=c(paste("Single, n=", n_s,sep=""), paste("NS, n=", n_g9)),
       col=c("red", "black"), bty = "n", lty=1, pch=pch_u)
dev.off()

pdf(paste(organ,"_deltaGC_exclude_CA_step=",step,"_d=",d,"-SI.pdf" ,sep=""), width=9.23, height = 6.52*1.5)
par(mfcol=c(3,1))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
# Brain multiple
organ="BN"
output = BN_output
g9_gc_e = output[[1]]
m_gc_e = output[[2]]
s_gc_e = output[[3]]
n_g9 = output[[4]] 
n_m = output[[5]] 
n_s = output[[6]]

plot(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", 
     #ylim=c(min(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]]), max(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]])),
     main=paste(organ, " High Allele - Low Allele, d=",g9_gc_e[[4]] ," bin size=",step,"_exclude maxTSN",sep=""),
     xlab="distance to maxTSN", ylab="delta GC", pch=pch_u, xlim=c(-d,d), las=1, frame=FALSE)

points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], g9_gc_e[[2]]-g9_gc_e[[3]], type="o", col="black", pch=pch_u)
points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=pch_u)
legend("topleft", legend=c(paste("Multiple, n=", n_m,sep=""), paste("NS, n=", n_g9)),
       col=c("blue", "black"), bty = "n", lty=1, pch=pch_u)



organ="LV"; 
output = LV_output
g9_gc_e = output[[1]]
m_gc_e = output[[2]]
s_gc_e = output[[3]]
n_g9 = output[[4]] 
n_m = output[[5]] 
n_s = output[[6]] 
### Single base plot
# width=9.23, height = 6.52)
#par(cex=1.5)
pch_u=15
#Liver Single base plot

plot(s_gc_e[[1]][1:length(s_gc_e[[2]])], s_gc_e[[2]]-s_gc_e[[3]], type="o", col="red", 
     #ylim=c(min(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]]), max(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]])),
     main=paste(organ, " High Allele - Low Allele, d=",g9_gc_e[[4]] ," bin size=",step,"_exclude maxTSN",sep=""),
     xlab="distance to maxTSN", ylab="delta GC", pch=pch_u, xlim=c(-d,d), las=1, frame=FALSE)
#abline(h=0)
#abline(v=0)



points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], g9_gc_e[[2]]-g9_gc_e[[3]], type="o", col="black", pch=pch_u)
#points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=19)

points(s_gc_e[[1]][1:length(s_gc_e[[2]])], s_gc_e[[2]]-s_gc_e[[3]], type="o", col="red", pch=pch_u)

legend("bottomleft", legend=c(paste("Single, n=", n_s,sep=""), paste("NS, n=", n_g9)),
       col=c("red", "black"), bty = "n", lty=1, pch=pch_u)

# Liver multiple base

plot(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", 
     #ylim=c(min(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]]), max(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]])),
     main=paste(organ, " High Allele - Low Allele, d=",g9_gc_e[[4]] ," bin size=",step,"_exclude maxTSN",sep=""),
     xlab="distance to maxTSN", ylab="delta GC", pch=pch_u, xlim=c(-d,d), las=1, frame=FALSE)

points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], g9_gc_e[[2]]-g9_gc_e[[3]], type="o", col="black", pch=pch_u)
points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=pch_u)
legend("topleft", legend=c(paste("Multiple, n=", n_m,sep=""), paste("NS, n=", n_g9)),
       col=c("blue", "black"), bty = "n", lty=1, pch=pch_u)

dev.off()

