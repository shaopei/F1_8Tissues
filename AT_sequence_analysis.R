setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/AT_sequence_analysis")

library("TmCalculator")
library(seqLogo)
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

SeqLogo <- function(seq, output, range=NULL) {
  #seq=m$HighAlleleSeq
  seq<- data.frame(do.call(rbind, strsplit(as.character(seq), "")))
  df <- seq_upperCase(seq)
  #create position weight matrix
  pwm <- apply(df, 1, proportion)
  if (!is.null(range)){
    p = makePWM((pwm[,range]))    
  }else{
    p = makePWM((pwm))  
  }
  
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


### check 2nd SNP, regardless what's in SNP1
organ="BN"
seq_a=read.table(paste(organ, "_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_SNP2_+-100_Long_ShortAlleleSeq.bed", sep=""))
dim(seq_a)
AT_1stbp_LongAllele=SeqLogo(seq_a$V7, paste(organ, "_AT_2ndSNP_LongAllele.pdf", sep=""))
AT_1stbp_ShortAllele=SeqLogo(seq_a$V8, paste(organ, "_AT_2ndSNP_ShortAllele.pdf", sep=""))

acgt_col=c("dark green", "blue", "orange" , "red")
acgt=c("A","C","G","T")

#pdf(paste(organ,"_deltaATCG_single.pdf" ,sep=""), width =7, height = 7)
#par(mfcol=c(2,1))
#par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
#par(mgp=c(3,1,0))
#par(cex.lab=2.2, cex.axis=2.2)
#par(cex=1.5)
pch_u=1
#d=s_gc[[4]]
d=100
w1=30
w2=30
s=0
range=(d+1-w1): (d+1+w2)
s_l=AT_1stbp_LongAllele
s_h=AT_1stbp_ShortAllele
plot((-w1+s):(w2+s),s_h[1,range] - s_l[1,range], col = acgt_col[1], type="o",
     ylim=c(-0.2,0.2), pch=pch_u,
     ylab="ShortAllele - LongAllele",
     xlab="Distance to maxTSN",
     main=organ,
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points((-w1+s):(w2+s),s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o", pch=pch_u)
}
legend("topleft", legend=acgt,
       col=acgt_col, bty = "n", lty=1, pch=pch_u)

### check 2nd SNP, SNP1 short allele ==C, long allele != C
organ="BN"
seq_a=read.table(paste(organ, "_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP2_+-100_Long_ShortAlleleSeq.bed", sep=""))
dim(seq_a)
AT_1stbp_LongAllele=SeqLogo(seq_a$V9, paste(organ, "_AT_2ndSNP_LongAllele.pdf", sep=""))
AT_1stbp_ShortAllele=SeqLogo(seq_a$V10, paste(organ, "_AT_2ndSNP_ShortAllele.pdf", sep=""))

acgt_col=c("dark green", "blue", "orange" , "red")
acgt=c("A","C","G","T")

#pdf(paste(organ,"_deltaATCG_single.pdf" ,sep=""), width =7, height = 7)
#par(mfcol=c(2,1))
#par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
#par(mgp=c(3,1,0))
#par(cex.lab=2.2, cex.axis=2.2)
#par(cex=1.5)
pch_u=1
#d=s_gc[[4]]
d=100
w1=30
w2=30
s=0
range=(d+1-w1): (d+1+w2)
s_l=AT_1stbp_LongAllele
s_h=AT_1stbp_ShortAllele
plot((-w1+s):(w2+s),s_h[1,range] - s_l[1,range], col = acgt_col[1], type="o",
     ylim=c(-0.20,0.20), pch=pch_u,
     ylab="ShortAllele - LongAllele",
     xlab="Distance to maxTSN",
     main=organ,
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points((-w1+s):(w2+s),s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o", pch=pch_u)
}
legend("topleft", legend=acgt,
       col=acgt_col, bty = "n", lty=1, pch=pch_u)


### check 1st SNP
organ="SP"
seq_a=read.table(paste(organ, "_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_+-30_Long_ShortAlleleSeq.bed", sep=""))
dim(seq_a)
AT_1stbp_LongAllele=SeqLogo(seq_a$V8, paste(organ, "_AT_1stbp_LongAllele.pdf", sep=""))
AT_1stbp_ShortAllele=SeqLogo(seq_a$V9, paste(organ, "_AT_1stbp_ShortAllele.pdf", sep=""))

acgt_col=c("dark green", "blue", "orange" , "red")
acgt=c("A","C","G","T")

#pdf(paste(organ,"_deltaATCG_single.pdf" ,sep=""), width =7, height = 7)
#par(mfcol=c(2,1))
#par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
#par(mgp=c(3,1,0))
#par(cex.lab=2.2, cex.axis=2.2)
#par(cex=1.5)
pch_u=1
#d=s_gc[[4]]
d=30
w1=0
w2=0
s=1
range=(d+1-w1): (d+1+w2)
s_l=AT_1stbp_LongAllele
s_h=AT_1stbp_ShortAllele
plot((-w1+s):(w2+s),s_h[1,range] - s_l[1,range], col = acgt_col[1], type="o",
     xlim=c(0,5),
     ylim=c(-0.12,0.12), pch=pch_u,
     ylab="ShortAllele - LongAllele",
     xlab="Distance to maxTSN",
     main=organ,
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points((-w1+s):(w2+s),s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o", pch=pch_u)
}
legend("topleft", legend=acgt,
       col=acgt_col, bty = "n", lty=1, pch=pch_u)

# check SNP2 when SNP1 is C in Short ATC in Long allele
organ="BN"
seq_a=read.table(paste(organ, "_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP2_+-100_Long_ShortAlleleSeq.bed", sep=""))
dim(seq_a)
AT_1stbp_LongAllele_SNP2=SeqLogo(seq_a$V9, paste(organ, "_AT_1stBpShortC_LongAllele_SNP2.pdf", sep=""))
AT_1stbp_ShortAllele_SNP2=SeqLogo(seq_a$V10, paste(organ, "_AT_1stBpShortC_ShortAllele_SNP2.pdf", sep=""))

d=100
w1=0
w2=0
s=2
range=(d+1-w1): (d+1+w2)
s_l=AT_1stbp_LongAllele_SNP2
s_h=AT_1stbp_ShortAllele_SNP2

for (i in 4:1){
  points((-w1+s):(w2+s),s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o", pch=pch_u)
}

# check SNP3 when SNP1 is C in Short ATC in Long allele
seq_a=read.table(paste(organ, "_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP3_+-100_Long_ShortAlleleSeq.bed", sep=""))
dim(seq_a)
AT_1stbp_LongAllele_SNP3=SeqLogo(seq_a$V9, paste(organ, "_AT_1stBpShortC_LongAllele_SNP3.pdf", sep=""))
AT_1stbp_ShortAllele_SNP3=SeqLogo(seq_a$V10, paste(organ, "_AT_1stBpShortC_ShortAllele_SNP3.pdf", sep=""))

s=3
s_l=AT_1stbp_LongAllele_SNP3
s_h=AT_1stbp_ShortAllele_SNP3

for (i in 4:1){
  points((-w1+s):(w2+s),s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o", pch=pch_u)
}

# check SNP4 when SNP1 is C in Short ATC in Long allele
seq_a=read.table(paste(organ, "_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP4_+-100_Long_ShortAlleleSeq.bed", sep=""))
dim(seq_a)
AT_1stbp_LongAllele_SNP4=SeqLogo(seq_a$V9, paste(organ, "_AT_1stBpShortC_LongAllele_SNP4.pdf", sep=""))
AT_1stbp_ShortAllele_SNP4=SeqLogo(seq_a$V10, paste(organ, "_AT_1stBpShortC_ShortAllele_SNP4.pdf", sep=""))

s=4
s_l=AT_1stbp_LongAllele_SNP4
s_h=AT_1stbp_ShortAllele_SNP4

for (i in 4:1){
  points((-w1+s):(w2+s),s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o", pch=pch_u)
}


# check SNP5 when SNP1 is C in Short ATC in Long allele
seq_a=read.table(paste(organ, "_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBpShortC_SNP5_+-100_Long_ShortAlleleSeq.bed", sep=""))
dim(seq_a)
AT_1stbp_LongAllele_SNP5=SeqLogo(seq_a$V9, paste(organ, "_AT_1stBpShortC_LongAllele_SNP5.pdf", sep=""))
AT_1stbp_ShortAllele_SNP5=SeqLogo(seq_a$V10, paste(organ, "_AT_1stBpShortC_ShortAllele_SNP5.pdf", sep=""))

s=5
s_l=AT_1stbp_LongAllele_SNP5
s_h=AT_1stbp_ShortAllele_SNP5

for (i in 4:1){
  points((-w1+s):(w2+s),s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o", pch=pch_u)
}

### use this one !!!
## check SNP1 to SNP5 ####################
organ="SP"
seq_a=read.table(paste(organ, "_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_+-30_Long_ShortAlleleSeq.bed", sep=""))
dim(seq_a)
AT_1stbp_LongAllele=SeqLogo(seq_a$V8, paste(organ,"_AT_1stbp_LongAllele.pdf",sep=""))
AT_1stbp_ShortAllele=SeqLogo(seq_a$V9, paste(organ,"_AT_1stbp_ShortAllele.pdf",sep=""))

acgt_col=c("dark green", "blue", "orange" , "red")
acgt=c("A","C","G","T")
acgt_SminusL <-NULL 
pch_u=1
d=30
w1=0
w2=0
s=1
range=(d+1-w1): (d+1+w2)
s_h=AT_1stbp_ShortAllele
s_l=AT_1stbp_LongAllele
plot((-w1+s):(w2+s),s_h[1,range] - s_l[1,range], col = acgt_col[1], type="o",
     xlim=c(0,5),
     ylim=c(-0.12,0.12), pch=pch_u,
     ylab="ShortAllele - LongAllele",
     xlab="SNP",
     main=organ,
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points((-w1+s):(w2+s),s_h[i,range] - s_l[i,range], col = acgt_col[i], type="o", pch=pch_u)
  acgt_SminusL[[acgt[i]]] <-list(s_h[i,range] - s_l[i,range])
  }
legend("topright", legend=acgt,
       col=acgt_col, bty = "n", lty=1, pch=pch_u)

d=100
w1=0
w2=0
range=(d+1-w1): (d+1+w2)
for (s in 2:5){
  #read.table("LV_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_SNP2_+-100_Long_ShortAlleleSeq.bed")
  seq_a=read.table(paste(organ, "_AT_4tunitIntersectNativeHMM_intersectRegion_strain_1stBp_SNP",s,"_+-100_Long_ShortAlleleSeq.bed", sep=""))
  cat ("SNP", s , dim(seq_a),"\n")
  LongAllele=SeqLogo(seq_a$V9, paste(organ, "_AT_1stBpShortC_LongAllele_SNP",s,".pdf", sep=""))
  ShortAllele=SeqLogo(seq_a$V10, paste(organ, "_AT_1stBpShortC_ShortAllele_SNP",s,".pdf", sep=""))
  shortMinusLong= ShortAllele - LongAllele
  
  for (i in 4:1){
    points((-w1+s):(w2+s), shortMinusLong[i,range] , col = acgt_col[i], type="o", pch=pch_u)
    acgt_SminusL[[acgt[i]]] <-list(unlist(acgt_SminusL[[acgt[i]]]), shortMinusLong[i,range])
      }
  
}
y=0.15
plot(1:5,unlist(acgt_SminusL[[acgt[1]]]), col = acgt_col[1], type="o",
     xlim=c(1,5),
     ylim=c(-y,y), pch=pch_u,
     ylab="ShortAllele - LongAllele",
     xlab="SNP",
     main=organ,
     las=1, frame=FALSE
)
abline(h=0, col="gray")
for (i in 4:1){
  points(1:5,unlist(acgt_SminusL[[acgt[i]]]), col = acgt_col[i], type="o", pch=pch_u)
}
legend("topright", legend=acgt,
       col=acgt_col, bty = "n", lty=1, pch=pch_u)


## combine A with T, C with G
#acgt=c("A","C","G","T")
plot(1:5,unlist(acgt_SminusL[[acgt[1]]])+unlist(acgt_SminusL[[acgt[4]]]), col = acgt_col[1], type="o",
     xlim=c(0,5),
     ylim=c(-0.15,0.15), pch=pch_u,
     ylab="ShortAllele - LongAllele",
     xlab="SNP",
     main=organ,
     las=1, frame=FALSE
)
abline(h=0, col="gray")
points(1:5,unlist(acgt_SminusL[[acgt[2]]])+unlist(acgt_SminusL[[acgt[3]]]), col = acgt_col[2], type="o", pch=pch_u)
legend("topleft", legend=c("AT", "CG"),
       col=c("dark green", "blue"), bty = "n", lty=1, pch=pch_u)





#######################

d=30
bin=13
seq_df=read.table("BN_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1_+-30_Long_ShortAlleleSeq.bed")
dim(seq_df)
for (bin in seq(8,13)){
range1=(d-bin*2+1):(d-bin)
range2=(d-bin+1):d
range3=(d+2):(d+bin+1)
range4=(d+bin+2):(d+bin*2+1)

for (i in 1:NROW(seq_df)){
  a=s2c(as.character(seq_df$V8[i]))
  a1=a[range1]
  a2=a[range2]
  a3=a[range3]
  a4=a[range4]
  seq_df$LGC_range1[i]= GC(a1)
  seq_df$LGC_range2[i]= GC(a2)
  seq_df$LGC_range3[i]= GC(a3)
  seq_df$LGC_range4[i]= GC(a4)
}
for (i in 1:NROW(seq_df)){
  a=s2c(as.character(seq_df$V9[i]))
    a1=a[range1]
  a2=a[range2]
  a3=a[range3]
  a4=a[range4]
    seq_df$SGC_range1[i]= GC(a1)
  seq_df$SGC_range2[i]= GC(a2)
  seq_df$SGC_range3[i]= GC(a3)
  seq_df$SGC_range4[i]= GC(a4)
}

vioplot(seq_df$SGC_range1, seq_df$SGC_range2, seq_df$SGC_range3, seq_df$SGC_range4,
        main= paste("short ","bin size = ", bin, sep=""),
        ylab="GC content", ylim=c(0,100), las=1, frame.plot = F)
}

vioplot(seq_df$LGC_range1, seq_df$LGC_range2, seq_df$LGC_range3, seq_df$LGC_range4,
        main= paste("long ","bin size = ", bin, sep=""),
        ylab="GC content", ylim=c(0,100), las=1, frame.plot = F)


vioplot(seq_df$LGC_range2, seq_df$SGC_range2, #seq_df$GC_range3, 
        main= paste("long short beforeAT","bin size = ", bin, sep=""),
        ylab="GC content", ylim=c(0,100), las=1, frame.plot = F)

vioplot(seq_df$LGC_range3, seq_df$SGC_range3, #seq_df$GC_range3, 
        main= paste("long short fromAT","bin size = ", bin, sep=""),
        ylab="GC content", ylim=c(0,100), las=1, frame.plot = F)


library("vioplot")
# panel57_GC_content_around_maxPause_noduplicate
# Figure4E, n=3456
dim(seq_df)
pdf("GC_content_around_maxPause_noduplicate.pdf", width=7, height = 3.5, useDingbats=FALSE)
par(mfrow=c(1,3))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)

vioplot(seq_df$GC_range2, seq_df$GC_range3, #seq_df$GC_range3, 
        main= paste("long ","bin size = ", bin, sep=""),
        ylab="GC content", ylim=c(0,100), las=1, frame.plot = F)




vioplot(seq_df$GC_range2, seq_df$GC_range3, 
        main= paste("short ","bin size = ", bin, sep=""),
        ylab="GC content", ylim=c(0,100), las=1, frame.plot = F)
abline(h=45, col="green")








name="BN"
name="LV"
SNPsAbundanceAround1npBed <-function(d=50, step=1,times=1, use.log=FALSE, use.sum=FALSE, name="", 
                                          OnlyAsTSS= FALSE, OnlynonAsTSSwithAsMaxTSN= FALSE, p_value_cut=NULL ,
                                          draw_plot=TRUE,
                                          s_l=NULL, m_l=NULL){ 
  d=10
  OneBp = read.table(paste(name,"_AT_4tunitIntersectNativeHMM_intersectRegion_strain_temp1.bed",sep=""), header = F)
  colnames(OneBp)[1:3]=c("chrm","chrmStart", "chrmEnd")

  g1 =  OneBp
  dim(g1); 
  
  g1$chrmStart = g1$chrmStart-d
  g1$chrmEnd = g1$chrmEnd+d

  # get the sum of SNPs in each bin (size = step)
  g1.SNPs <- read_read_mat_SNPs (SNP.bw, g1[,c(1:6)], step, times=times, use.log=use.log)

 show.window <- d
  if (use.sum){
    a = colSums(g1.SNPs)
  }else{
    a = colMeans(g1.SNPs)
  }
  
  if (step >1){
    x <- seq(-1*(show.window)+step/2, show.window, step)
  }else{
    x <- seq(-1*(show.window), show.window, step)
  }
  
  if(draw_plot){
    #pdf(paste(name,"_AsTSS_SNPs_d=",d,"_step=",step,"-S.pdf",sep=""), width=9.23, height = 6.52)
    par(mfcol=c(2,1))
    par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
    par(mgp=c(3,1,0))
    par(cex.lab=2.2, cex.axis=2.2)
    #par(cex=1.5)
    pch_u=15
    plot(x,a , type="o") #,ylim=c(0,0.2)
    points(x,a, type="o", col="red")
    
    plot(x, a,col="black",  xlab="distance to maxTSN", ylab="SNPs mean", 
         main=name, type="o", ylim=c(0,max(a)), pch=pch_u, las=1, frame=FALSE)
    points(x, s_plot, col="red", type="o", pch=pch_u)
    
    if (! is.null(s_l)){
      for (l in s_l){
        abline(v=l-5, col="green")
        abline(v=l+4, col="orange")
      }
    }
    
    legend("topleft", legend=c(paste("Single, n=", dim(s)[1],sep=""), paste("NS, n=", dim(g9)[1])),
           col=c("red", "black"), bty = "n", lty=1, pch=pch_u)
    
    plot(x, m_plot, col="blue", xlab="distance to maxTSN", ylab="SNPs mean", main="", type="o", 
         ylim=c(0,max(a,b,s_plot,m_plot)), las=1, pch=pch_u, frame=FALSE)
    points(x, b,col="black", type="o", pch=pch_u)
    
    if (! is.null(m_l)){
      for (l in m_l){
        abline(v=l-5, col="green")
        abline(v=l+4, col="orange")
      }
    }
    
    
    
    legend("topleft", legend=c(paste("Multiple, n=", dim(m)[1],sep=""), paste("NS, n=", dim(g9)[1])),
           col=c("blue", "black"), bty = "n", lty=1, pch=pch_u) #, cex=2)
    
    dev.off()
  }
  return (list(x,s_plot, m_plot, dim(s), dim(m), b, dim(g9)))
  
}
