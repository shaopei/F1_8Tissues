
library("TmCalculator")
library(seqLogo)
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/GC_content")


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
  
  
  g9_gc_e=High_Low_GC_excludeCA(g9, name=paste(organ,"_fdr>0.9_d=",d,sep=""), step=step)
  m_gc_e=High_Low_GC_excludeCA(m, name=paste(organ,"_m_d=",d,sep=""), step=step)
  s_gc_e=High_Low_GC_excludeCA(s, name=paste(organ,"_s_d=",d,sep=""), step=step)
  return(list(g9_gc_e, m_gc_e, s_gc_e, dim(g9)[1], dim(m)[1], dim(s)[1] ))
}

d=35; step=5
organ="LV"; 
LV_df=read.table(paste(organ, "_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs_+-",d,"_High_LowAlleleSeq.bed", sep=""))
LV_output = df_process(LV_df, organ, step)

organ="BN"
BN_df=read.table(paste(organ, "_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs_+-",d,"_High_LowAlleleSeq.bed", sep=""))
#name=organ

BN_output = df_process(BN_df, organ, step)


### Single base plot
g9_gc_e = LV_output[[1]]
m_gc_e = LV_output[[2]]
s_gc_e = LV_output[[3]]

#pdf(paste(organ,"_deltaGC_exclude_CA_step=",step,"_d=",d,".pdf" ,sep=""), width = 10, height = 10)
plot(s_gc_e[[1]][1:length(s_gc_e[[2]])], s_gc_e[[2]]-s_gc_e[[3]], type="o", col="pink", 
     ylim=c(min(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]]), max(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]])),
     main=paste("High Allele - Low Allele, d=",g9_gc_e[[4]] ," bin size=",step,"_exclude maxTSN",sep=""),
     xlab="distance to maxTSN", ylab="delta GC", pch=6, xlim=c(-d,d))
#abline(h=0)
#abline(v=0)


points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], g9_gc_e[[2]]-g9_gc_e[[3]], type="o", col="gray", pch=6)
#points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=19)

points(s_gc_e[[1]][1:length(s_gc_e[[2]])], s_gc_e[[2]]-s_gc_e[[3]], type="o", col="pink", pch=6)




g9_gc_e = BN_output[[1]]
m_gc_e = BN_output[[2]]
s_gc_e = BN_output[[3]]



points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], g9_gc_e[[2]]-g9_gc_e[[3]], type="o", col="gray", pch=1)
#points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=19)

points(s_gc_e[[1]][1:length(s_gc_e[[2]])], s_gc_e[[2]]-s_gc_e[[3]], type="o", col="pink", pch=1)


a_g9_gc_e = ((BN_output[[1]][[2]]- BN_output[[1]][[3]]) + (LV_output[[1]][[2]]- LV_output[[1]][[3]]))/2
a_m_gc_e =((BN_output[[2]][[2]]- BN_output[[2]][[3]]) + (LV_output[[2]][[2]]- LV_output[[2]][[3]]))/2
a_s_gc_e = ((BN_output[[3]][[2]]- BN_output[[3]][[3]]) + (LV_output[[3]][[2]]- LV_output[[3]][[3]]))/2


#points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], a_g9_gc_e, type="o", col="black", pch=19)
plot(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], a_g9_gc_e, type="o", col="black", pch=15,
     ylim=c(min(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]]), max(s_gc_e[[2]]-s_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]])),
     main=paste("High Allele - Low Allele, d=",g9_gc_e[[4]] ," bin size=",step,"_exclude maxTSN",sep=""),
     xlab="distance to maxTSN", ylab="delta GC",
     xlim=c(-d,d))
#points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=19)

points(s_gc_e[[1]][1:length(s_gc_e[[2]])], a_s_gc_e, type="o", col="red", pch=15)

legend("bottomleft", legend=c("Single", "NS"),
       col=c("red", "black"), bty = "n", lty=1, pch=c(15,15))

legend("bottomleft", legend=c("Single Mean", "Single Brain", "Single Liver", "NS Mean", "NS Brain", "NS Liver"),
       col=c("red","pink","pink", "black","gray","gray"), bty = "n", lty=1, pch=c(19,1,6,19,1,6))

dev.off()


# fisher's exact test
#BN_output 
test_output = LV_output #  return(list(g9_gc_e, m_gc_e, s_gc_e, dim(g9)[1], dim(m)[1], dim(s)[1] ))
g9_gc_e = test_output[[1]] #  return (list(x, HighTm, LowTm, show.window))
m_gc_e = test_output[[2]]
s_gc_e = test_output[[3]]
g9_len=test_output[[4]]
m_len=test_output[[5]]
s_len=test_output[[6]]

p.value <- NULL
odds.ratio <- NULL
testors <- NULL
for (i in (1:length(g9_gc_e[[2]]))){
# testor = rbind(c(round(s_gc_e[[2]][i]*s_len*step/100), round(g9_gc_e[[2]][i]*g9_len*step/100)),
#   c(round(s_gc_e[[3]][i]*s_len*step/100),round(g9_gc_e[[3]][i]*g9_len*step/100) )
#                ); testor
testor <-  matrix(c(round(g9_gc_e[[2]][i]*g9_len*step/100), round(g9_gc_e[[3]][i]*g9_len*step/100),
                    round(s_gc_e[[2]][i]*s_len*step/100), round(s_gc_e[[3]][i]*s_len*step/100)),
         nrow = 2,
         dimnames = list(Allele = c("High", "Low"),
                         TSS = c("NS", "Single"))); testor
# multiple
if(0){
testor <-  matrix(c(round(g9_gc_e[[2]][i]*g9_len*step/100), round(g9_gc_e[[3]][i]*g9_len*step/100),
                    round(m_gc_e[[2]][i]*m_len*step/100), round(m_gc_e[[3]][i]*m_len*step/100)),
                  nrow = 2,
                  dimnames = list(Allele = c("High", "Low"),
                                  TSS = c("NS", "Multiple"))); testor
}
testors=rbind(testors, as.vector(t(testor)))
#cat(testor, "\n")
f = fisher.test(testor); f
p.value= c(p.value, f$p.value)
odds.ratio = c(odds.ratio, f$estimate)
}

plot(g9_gc_e[[1]],p.value, type="o" , xlab="dist to maxTSN")

# Multiple base plot



g9_gc_e = BN_output[[1]]
m_gc_e = BN_output[[2]]
s_gc_e = BN_output[[3]]

#pdf(paste(organ,"_deltaGC_exclude_CA_step=",step,"_d=",d,".pdf" ,sep=""), width = 10, height = 10)
plot(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="light blue", 
     #ylim=c(min(m_gc_e[[2]]-m_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]]), max(m_gc_e[[2]]-m_gc_e[[3]], m_gc_e[[2]]-m_gc_e[[3]])),
     ylim=c(-0.9, 0.5),
     main=paste("High Allele - Low Allele, d=",g9_gc_e[[4]] ," bin size=",step,"_exclude maxTSN",sep=""),
     xlab="distance to maxTSN", ylab="delta GC", pch=1,
     xlim=c(-d,d))

points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], g9_gc_e[[2]]-g9_gc_e[[3]], type="o", col="gray", pch=1)
#points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=19)

points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="light blue", pch=1)


g9_gc_e = LV_output[[1]]
m_gc_e = LV_output[[2]]
s_gc_e = LV_output[[3]]


#abline(h=0)
#abline(v=0)


points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], g9_gc_e[[2]]-g9_gc_e[[3]], type="o", col="gray", pch=6)
#points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=19)

points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="light blue", pch=6)




a_g9_gc_e = ((BN_output[[1]][[2]]- BN_output[[1]][[3]]) + (LV_output[[1]][[2]]- LV_output[[1]][[3]]))/2
a_m_gc_e =((BN_output[[2]][[2]]- BN_output[[2]][[3]]) + (LV_output[[2]][[2]]- LV_output[[2]][[3]]))/2
a_s_gc_e = ((BN_output[[3]][[2]]- BN_output[[3]][[3]]) + (LV_output[[3]][[2]]- LV_output[[3]][[3]]))/2

points(g9_gc_e[[1]][1:length(g9_gc_e[[2]])], a_g9_gc_e, type="o", col="black", pch=19)
#points(m_gc_e[[1]][1:length(m_gc_e[[2]])], m_gc_e[[2]]-m_gc_e[[3]], type="o", col="blue", pch=19)

points(m_gc_e[[1]][1:length(m_gc_e[[2]])], a_m_gc_e, type="o", col="blue", pch=19)

legend("bottomleft", legend=c("Multiple Mean", "Multiple Brain", "Multiple Liver", "NS Mean", "NS Brain", "NS Liver"),
       col=c("blue","light blue","light blue", "black","gray","gray"), bty = "n", lty=1, pch=c(19,1,6,19,1,6))

dev.off()


















### calculate the average of the delta GC at each base of each organ
Organ=c("BN","LV")#"HT", "SP", "GI", "ST")
g9_deltaGC <- NULL
s_deltaGC <- NULL
m_deltaGC <- NULL
d=35; step=1

for(t in Organ){
  organ=t; 
  df=read.table(paste(organ, "_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs_+-",d,"_High_LowAlleleSeq.bed", sep=""))
  output = df_process(df, organ, step)
  g9_gc_e = output[[1]]
  m_gc_e = output[[2]]
  s_gc_e = output[[3]]  # (list(x, HighTm, LowTm, show.window))
  
  g9_deltaGC = rbind(g9_deltaGC, (g9_gc_e[[2]] - g9_gc_e[[3]]))
  m_deltaGC = rbind(m_deltaGC, (m_gc_e[[2]] - m_gc_e[[3]]))
  s_deltaGC = rbind(s_deltaGC,(s_gc_e[[2]] - s_gc_e[[3]]))
    }

save.image("TSS_GC_NS.Single.Multiple.RData")
load("TSS_GC_NS.Single.Multiple.RData")

colMedian <- function(x){
  return(apply(x, 2, median))
}

step=5
sumByBin <- function(x, bin=step){
  # the last bin will be removed
  output_bin<- NULL
  for (i in (seq(1,length(x)-bin,bin))){
    output_bin = c(output_bin, sum(x[i:(i+bin-1)]))
  }
  return (output_bin)
}

meanByBin <- function(x, bin=step){
  # the last bin will be removed
  output_bin<- NULL
  for (i in (seq(1,length(x)-bin,bin))){
    output_bin = c(output_bin, mean(x[i:(i+bin-1)]))
  }
  return (output_bin)
}

g9_deltaGC_mean = sumByBin(colMeans(g9_deltaGC))
s_deltaGC_mean = sumByBin(colMeans(s_deltaGC))
m_deltaGC_mean = sumByBin(colMeans(m_deltaGC))

g9_deltaGC_mean = meanByBin(colMeans(g9_deltaGC))
s_deltaGC_mean = meanByBin(colMeans(s_deltaGC))
m_deltaGC_mean = meanByBin(colMeans(m_deltaGC))

show.window = floor(length(g9_gc_e[[1]])/2)
if (step >1){
  x <- seq(-1*(show.window)+floor(step/2), show.window, step)
}else{
  x <- seq(-1*(show.window), show.window, step)
}
x=x[1:length(g9_deltaGC_mean)]

plot(x,g9_deltaGC_mean, type="o", col="black", xlab="distance to maxTSN", ylab="sum(delta GC in bin)", pch=19,
       main=paste("High Allele - Low Allele, d=",d ," bin size=",step,"_exclude maxTSN",sep=""),
       xlim=c(-d,d), ylim=c(-2,1))
points(x,m_deltaGC_mean, type="o", col="blue", pch=19)
points(x, s_deltaGC_mean, type="o", col="red", pch=19)

abline(v=0, col="gray")
abline(v=-1, col="gray")



