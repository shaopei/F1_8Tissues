
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/GC_content")
organ="LV"; d=35; name=organ
df=read.table(paste(organ, "_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_BinomialTest_maxTSNs_+-",d,"_High_LowAlleleSeq.bed", sep=""))

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

write.table(unique(m[,1:6]), file=paste(organ, "_allReads_TSS_5mat5pat_uniq_AsTSSMultipleBaseDriven.bed", sep=""), 
            quote=F, sep="\t", row.names = F, col.names = F)
write.table(unique(s[,1:6]), file=paste(organ, "_allReads_TSS_5mat5pat_uniq_AsTSSSingleBaseDriven.bed", sep=""), 
            quote=F, sep="\t", row.names = F, col.names = F)


AT2GC_SNP <- function(df, name="" ,step=5, excludeCA=FALSE){
  AT2GC_SNP=c("AG","AC","TG","TC")
  TSS_with_AT2GC_SNP <- NULL
  TSS_with_other_SNP <- NULL
  TSS_without_SNP <- NULL
  for (i in 1:NROW(df)){
    i_TSS_with_AT2GC_SNP <- NULL
    i_TSS_with_other_SNP <- NULL
    i_TSS_without_SNP <- NULL
    high_a=s2c(as.character(df$HighAlleleSeq[i]))
    low_a=s2c(as.character(df$LowAlleleSeq[i]))
    len_a = length(high_a)
    d = (len_a-1)/2
    
    if (excludeCA){
      #CA at a[d:(d+1)]
      high_a = c(high_a[1:(d-1)] ,"N","N",high_a[(d+2):len_a])
      low_a=c(low_a[1:(d-1)] ,"N","N",low_a[(d+2):len_a])
    }
    for (j in seq(1,(len_a-step), step)){
      snp_a=FALSE
      snp_o=FALSE
      snp_n=FALSE
      
      sub_h=high_a[j:(j+step-1)] ; sub_l=low_a[j:(j+step-1)]
      # if there is seq difference between high and low allele
      if(sum(sub_h != sub_l)){  
        #cat (i, j, "\n")
        snp = which(sub_h != sub_l)
        for (k in 1:length(snp)){
          h2l_snp = c2s(c(sub_h[snp[k]], sub_l[snp[k]]))
          if (h2l_snp %in% AT2GC_SNP){
            # contain AT to GC snps
            snp_a=TRUE
          }
        }
        if (! snp_a){
          # other snp
          snp_o = TRUE
        }
      }else{
        # no snp
        snp_n=TRUE 
      }
      
      i_TSS_with_AT2GC_SNP <- c(i_TSS_with_AT2GC_SNP, snp_a)
      i_TSS_with_other_SNP <- c(i_TSS_with_other_SNP, snp_o)
      i_TSS_without_SNP <- c(i_TSS_without_SNP, snp_n)
    }
    TSS_with_AT2GC_SNP <- rbind(i_TSS_with_AT2GC_SNP, TSS_with_AT2GC_SNP)
    TSS_with_other_SNP <- rbind(i_TSS_with_other_SNP, TSS_with_other_SNP)
    TSS_without_SNP <- rbind(i_TSS_without_SNP, TSS_without_SNP)
  }
  show.window = floor(len_a/2)
  if (step >1){
    x <- seq(-1*(show.window)+floor(step/2), show.window, step)
  }else{
    x <- seq(-1*(show.window), show.window, step)
  }
  x=x[1:length(colSums(TSS_with_other_SNP))]
  plot(x, colSums(TSS_with_other_SNP), col="black", type="o", 
       xlab="distance to maxTSN", main=name,
       ylim = c(0,max(colSums(TSS_with_other_SNP),colSums(TSS_with_AT2GC_SNP))) )   
  points(x, colSums(TSS_with_AT2GC_SNP), col="red", type="o")   
  legend("topleft", legend=c("TSS_with_AT2GC_SNP", "TSS_with_other_SNP"),
         col=c("red", "black"), bty = "n", lty=1, pch=1)
  
  return (list(x,colSums(TSS_with_AT2GC_SNP), colSums(TSS_with_other_SNP), colSums(TSS_without_SNP)))
}

#dev.off()
par(mfrow=c(3,1))
step=5
#g9_tss = AT2GC_SNP(g9, step=step, "fdr>0.9, include CA")
#s_tss = AT2GC_SNP(s, step=step, "single, include CA")
#m_tss=AT2GC_SNP(m, step=step, "multiple, inlcude CA")


g9_tss = AT2GC_SNP(g9, step=step, "fdr>0.9, exclude CA", excludeCA=TRUE)
s_tss = AT2GC_SNP(s, step=step, "single, exclude CA", excludeCA=TRUE)
m_tss = AT2GC_SNP(m, step=step, "multiple, exlcude CA", excludeCA=TRUE)

if (0){
  fdr_cutoff=0.05
  p.value <- NULL
  odds.ratio <- NULL
  testors <- NULL
  
  for (n in (1:length(s_tss[[2]]))){
    
    testor <-  matrix(c(s_tss[[2]][n], s_tss[[4]][n],
                        g9_tss[[2]][n], g9_tss[[4]][n]),
                      nrow = 2,
                      dimnames = list(TSS = c("With ATtoGC SNPs", "With other SNPs"),
                                      KS.Test = c("Single", "NS"))); testor
    
    f = fisher.test(testor, alternative = "greater" ); f
    p.value= c(p.value, f$p.value)
    odds.ratio = c(odds.ratio, f$estimate)
  }
  adjust.p = p.adjust(p.value, method="fdr")
  par(mfrow=c(4,1))
  #plot(s_tss[[1]],p.value, type="o" , xlab="dist to maxTSN", main=paste(organ,"_SingleBaseTSS_SNPs_bin=",step, sep=""), las=1, col="red")
  plot(s_tss[[1]], -log10(adjust.p), type="o" , xlab="dist to maxTSN", main=" TSS_with_AT2GC_SNP / TSS_without_SNP", las=1, col="red")
  abline(h=-1*log10(fdr_cutoff),col="gray")
  plot(s_tss[[1]],odds.ratio, type="o" , xlab="dist to maxTSN",  main=paste(organ,"_SingleBaseTSS_SNPs_bin=",step, sep=""), las=1)
  abline(h=1,col="gray")
  # plot(s_tss[[1]], p.value, type="o" , 
  #      ylab="unadjust p-value",
  #      xlab="dist to maxTSN",main=" TSS_with_AT2GC_SNP / TSS_without_SNP", las=1, col="red")
  # abline(h=0.05,col="gray")
  
  
  
  p.value <- NULL
  odds.ratio <- NULL
  testors <- NULL
  
  for (n in (1:length(s_tss[[2]]))){
    
    testor <-  matrix(c(s_tss[[3]][n], s_tss[[4]][n],
                        g9_tss[[3]][n], g9_tss[[4]][n]),
                      nrow = 2,
                      dimnames = list(TSS = c("With ATtoGC SNPs", "With other SNPs"),
                                      KS.Test = c("Single", "NS"))); testor
    
    f = fisher.test(testor, alternative = "greater" ); f
    p.value= c(p.value, f$p.value)
    odds.ratio = c(odds.ratio, f$estimate)
  }
  adjust.p = p.adjust(p.value, method="fdr")
  
  plot(s_tss[[1]], -log10(adjust.p), type="o" , xlab="dist to maxTSN", 
       main=paste(organ,"_SingleBaseTSS_SNPs_bin=",step," TSS_with_other_SNP/ TSS_without_SNP", sep=""), las=1, col="red")
  abline(h=-1*log10(fdr_cutoff),col="gray")
  plot(s_tss[[1]],odds.ratio, type="o" , xlab="dist to maxTSN", main=paste(organ,"_SingleBaseTSS_SNPs", sep=""), las=1)
  abline(h=1,col="gray")
  
  
  
  
  #dev.off()
  par(mfrow=c(3,1))
  p.value <- NULL
  odds.ratio <- NULL
  testors <- NULL
  
  for (n in (1:length(s_tss[[2]]))){
    
    testor <-  matrix(c(s_tss[[2]][n], s_tss[[3]][n],
                        g9_tss[[2]][n], g9_tss[[3]][n]),
                      nrow = 2,
                      dimnames = list(TSS = c("With ATtoGC SNPs", "With other SNPs"),
                                      KS.Test = c("Single", "NS"))); testor
    
    f = fisher.test(testor, alternative = "greater" ); f
    p.value= c(p.value, f$p.value)
    odds.ratio = c(odds.ratio, f$estimate)
  }
  adjust.p = p.adjust(p.value, method="fdr")
  plot(s_tss[[1]], -log10(adjust.p), type="o" , xlab="dist to maxTSN", main=" TSS_with_AT2GC_SNP/ TSS_with_other_SNP", las=1, col="red")
  abline(h=-1*log10(fdr_cutoff),col="gray")
  plot(s_tss[[1]],odds.ratio, type="o" , xlab="dist to maxTSN", main=paste(organ,"_SingleBaseTSS_SNPs", sep=""), las=1)
  abline(h=1,col="gray")
  plot(s_tss[[1]], p.value, type="o" , 
       ylab="unadjust p-value",
       xlab="dist to maxTSN", main=" TSS_with_AT2GC_SNP/ TSS_with_other_SNP", las=1, col="red")
  abline(h=0.05,col="gray")
  text(s_tss[[1]][which(p.value <= 0.1)],p.value[which(p.value <= 0.1)], label=paste(round(p.value[which(p.value <= 0.1)], digits = 2), sep=" "))
  
}

par(mfrow=c(3,1))
p.value <- NULL
odds.ratio <- NULL
testors <- NULL

for (n in (1:length(s_tss[[2]]))){
  
  testor <-  matrix(c(s_tss[[2]][n], dim(s)[1],
                      g9_tss[[2]][n], dim(g9)[1]),
                    nrow = 2,
                    dimnames = list(TSS = c("TSS with ATtoGC SNPs", "All TSS"),
                                    KS.Test = c("Single", "NS"))); testor
  
  f = fisher.test(testor, alternative = "greater" ); f
  p.value= c(p.value, f$p.value)
  odds.ratio = c(odds.ratio, f$estimate)
}
adjust.p = p.adjust(p.value, method="fdr")
plot(s_tss[[1]], -log10(adjust.p), type="o" , xlab="dist to maxTSN", main=" TSS_with_AT2GC_SNP/ Total TSS", las=1, col="red")
abline(h=-1*log10(fdr_cutoff),col="gray")
plot(s_tss[[1]],odds.ratio, type="o" , xlab="dist to maxTSN", main=paste(organ,"_SingleBaseTSS_SNPs", sep=""), las=1)
abline(h=1,col="gray")
plot(s_tss[[1]], p.value, type="o" , 
     ylab="unadjust p-value",
     xlab="dist to maxTSN",  main=" TSS_with_AT2GC_SNP/ Total TSS", las=1, col="red")
abline(h=0.05,col="gray")
text(s_tss[[1]][which(p.value <= 0.1)],p.value[which(p.value <= 0.1)], label=paste(round(p.value[which(p.value <= 0.1)], digits = 2), sep=" "))


## multiple
if(0){
par(mfrow=c(3,1))
p.value <- NULL
odds.ratio <- NULL
testors <- NULL

for (n in (1:length(m_tss[[2]]))){
  
  testor <-  matrix(c(m_tss[[2]][n], dim(s)[1],
                      g9_tss[[2]][n], dim(g9)[1]),
                    nrow = 2,
                    dimnames = list(TSS = c("With ATtoGC SNPs", "With other SNPs"),
                                    KS.Test = c("Single", "NS"))); testor
  
  f = fisher.test(testor, alternative = "greater" ); f
  p.value= c(p.value, f$p.value)
  odds.ratio = c(odds.ratio, f$estimate)
}
adjust.p = p.adjust(p.value, method="fdr")
plot(m_tss[[1]], -log10(adjust.p), type="o" , xlab="dist to maxTSN", main=" TSS_with_AT2GC_SNP/ Total TSS", las=1, col="red")
abline(h=-1*log10(fdr_cutoff),col="gray")
plot(m_tss[[1]],odds.ratio, type="o" , xlab="dist to maxTSN", main=paste(organ,"_SingleBaseTSS_SNPs", sep=""), las=1)
abline(h=1,col="gray")
plot(m_tss[[1]], p.value, type="o" , 
     ylab="unadjust p-value",
     xlab="dist to maxTSN",  main=" TSS_with_AT2GC_SNP/ Total TSS", las=1, col="red")
abline(h=0.05,col="gray")
text(m_tss[[1]][which(p.value <= 0.1)],p.value[which(p.value <= 0.1)], label=paste(round(p.value[which(p.value <= 0.1)], digits = 2), sep=" "))
}


