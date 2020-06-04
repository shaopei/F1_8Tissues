# Rscript getSNPsAbundance.R

require(bigWig)


### SNP locations
#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation/SNPs_distribution")
#source("/Users/shaopei/Box\ Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/heatmap/heatmaps.R")

#SNP.bw <- "/Volumes/SPC_SD/KD_IGV/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"
SNP.bw <- "P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"

read_read_mat_SNPs <-function (SNP.bw , bed6, step=2, navg = 20, times=1, use.log=FALSE)
{
  bw.plus  <- load.bigWig(SNP.bw )
  
  hCountMatrix <- bed.step.bpQuery.bigWig(bw.plus, bed6[,c(1:3)] , step=step, abs.value=TRUE, op = "sum")
  hCountMatrix <- lapply(1:NROW(hCountMatrix), function(i){ if(bed6[i,6]=="-") return(rev(hCountMatrix[[i]])) else return(hCountMatrix[[i]])} );
  if (!use.log){
    hmat <- times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) ;
  } else {
    hmat <- log(times * matrix(unlist(hCountMatrix), nrow= NROW(bed6), byrow=TRUE) + 1) ;
  }
  #avgMat <- t(sapply(1:floor(NROW(hmat)/navg), function(x) {colMeans(hmat[((x-1)*navg+1):min(NROW(hmat),(x*navg)),])}))
  
  unload.bigWig(bw.plus);
  
  return(hmat);    
}

SNPsAbundanceAroundMaxTSN <-function(maxTSN, d=20, step=1,times=1, use.log=FALSE, use.sum=FALSE, name=""){ 
#d=20
#maxTSN=read.table(file = "BN_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr1.bed", header = T)
g1 = maxTSN[maxTSN$fdr<=0.1,c(1:4,11,10)]
g9 = maxTSN[maxTSN$fdr>=0.9,c(1:4,11,10)]

dim(g1); 
dim(g9)

g1$chrmStart = g1$chrmStart-d
g1$chrmEnd = g1$chrmEnd+d
g9$chrmStart = g9$chrmStart-d
g9$chrmEnd = g9$chrmEnd+d

  # get the sum of SNPs in each bin (size = step)
  g1.SNPs <- read_read_mat_SNPs (SNP.bw, g1[,c(1:6)], step, times=times, use.log=use.log)
  g9.SNPs <- read_read_mat_SNPs (SNP.bw, g9[,c(1:6)], step, times=times, use.log=use.log)
  
show.window <- d


  if (use.sum){
    a = colSums(g1.SNPs)
    b = colSums(g9.SNPs)
  }else{
    a = colMeans(g1.SNPs)
    b = colMeans(g9.SNPs)
  }
  
  if (step >1){
    x <- seq(-1*(show.window)+step/2, show.window, step)
  }else{
    x <- seq(-1*(show.window), show.window, step)
  }
  
  pdf(paste(name,"SNPs.pdf",sep="_"), width=6, height = 6)
  par(mfrow=c(2,1))
  #par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
  #par(mgp=c(3,1,0))
  #par(cex.lab=2.2, cex.axis=2.2)  
  
  plot(x, b, col="blue", xlab="distance to maxTSN", ylab="SNPs mean", main=name, type="o", ylim=c(0,max(a,b)), las=1)
  points(x, a,col="red", type="o")
  legend("topright", legend=c("fdr<=0.1", "fdr>0.9"),
         col=c("red", "blue"), bty = "n", lty=1, pch=1)
  plot(x,a-b, type="o", xlab="distance to maxTSN",ylab="substract", las=1, main="substract", ylim=c(-0.05,0.15))
  dev.off()

}

for(t in c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")){
df=read.table(file = paste(t,"_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr1.bed",sep = ""), header = T)
SNPsAbundanceAroundMaxTSN(df, name=t)
}


SNPsAbundanceAroundMaxTSNInTSS <-function(d=50, step=1,times=1, use.log=FALSE, use.sum=FALSE, name="", OnlyAsTSS= FALSE, OnlynonAsTSSwithAsMaxTSN= FALSE, p_value_cut=NULL , draw_plot=TRUE){ 
  #d=202
  AllTSS_maxTSN = read.table(paste(name,"_allReads_TSS_5mat5pat_uniq_pValue_maxTSNs.bed",sep=""), header = F)
  colnames(AllTSS_maxTSN)[7:10]=c("masked_p_value","masked_fdr", "unmasked_p_value","unmasked_fdr")
  colnames(AllTSS_maxTSN)[11:13]=c("chrm","chrmStart", "chrmEnd")
  # TSS not allelic different (fdr>0.9)
  g1 =  AllTSS_maxTSN[AllTSS_maxTSN$unmasked_fdr <= 0.1,c(11:16,9)]
  g9=  AllTSS_maxTSN[AllTSS_maxTSN$unmasked_fdr > 0.9,c(11:16,9)]
  dim(g1); max(g1$unmasked_p_value); min(g1$unmasked_p_value)
  dim(g9); max(g9$unmasked_p_value); min(g9$unmasked_p_value)
  
  g1$chrmStart = g1$chrmStart-d
  g1$chrmEnd = g1$chrmEnd+d
  g9$chrmStart = g9$chrmStart-d
  g9$chrmEnd = g9$chrmEnd+d
  
  # get the sum of SNPs in each bin (size = step)
  g1.SNPs <- read_read_mat_SNPs (SNP.bw, g1[,c(1:6)], step, times=times, use.log=use.log)
  g9.SNPs <- read_read_mat_SNPs (SNP.bw, g9[,c(1:6)], step, times=times, use.log=use.log)
  
  #AsTSS_maxTSN=read.table(file = paste(name,"_allReads_TSS_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue_maxTSNs.bed",sep=""), header = F)
  AsTSS_maxTSN=read.table(file = paste(name,"_allReads_TSS_5mat5pat_uniq_maskedVSunmasked_pValue_maxTSNs.bed",sep=""), header = F)
  
  colnames(AsTSS_maxTSN)[7:10]=c("masked_p_value","masked_fdr", "unmasked_p_value","unmasked_fdr")
  colnames(AsTSS_maxTSN)[11:13]=c("chrm","chrmStart", "chrmEnd")
  
  if (OnlyAsTSS){  # only use as.TSS. Remove TSS with as.maxTSN
    AsTSS_maxTSN=AsTSS_maxTSN[AsTSS_maxTSN$unmasked_fdr <= 0.1 ,]
  } else if (OnlynonAsTSSwithAsMaxTSN){  # only nonAS.TSS with as.maxTSN
    AsTSS_maxTSN=AsTSS_maxTSN[AsTSS_maxTSN$unmasked_fdr > 0.1 ,]
  }
  if (is.null(p_value_cut)){
    p_value_cut = max(AsTSS_maxTSN$unmasked_p_value[AsTSS_maxTSN$unmasked_fdr<=0.1])
  }
  
  # asTSS driven by multiple base
  m = AsTSS_maxTSN[AsTSS_maxTSN$masked_p_value <= p_value_cut,c(11:16,7)]
  # asTSS driven by single base
  s = AsTSS_maxTSN[AsTSS_maxTSN$masked_p_value > p_value_cut,c(11:16,7)]
  
  # if (dim(s)[1] > dim(m)[1]){
  # s = s[order(s$masked_p_value, decreasing = TRUE),][1:dim(m)[1],c(11:16,7)]
  # }else{
  #   s = s[,c(11:16,7)]
  # }
  
  cat(name, "\n")
  cat ("dim(m), max(m$masked_p_value), min(m$masked_p_value)", "\n")
  cat (dim(m), max(m$masked_p_value), min(m$masked_p_value), "\n")
  cat("  dim(s); max(s$masked_p_value); min(s$masked_p_value)", "\n")
  cat(dim(s), max(s$masked_p_value), min(s$masked_p_value), "\n")
  
  
  s$chrmStart = s$chrmStart-d
  s$chrmEnd = s$chrmEnd+d
  m$chrmStart = m$chrmStart-d
  m$chrmEnd = m$chrmEnd+d
  
  # get the sum of SNPs in each bin (size = step)
  s.SNPs <- read_read_mat_SNPs (SNP.bw, s[,c(1:6)], step, times=times, use.log=use.log)
  m.SNPs <- read_read_mat_SNPs (SNP.bw, m[,c(1:6)], step, times=times, use.log=use.log)
  
  
  show.window <- d
  if (use.sum){
    a = colSums(g1.SNPs)
    b = colSums(g9.SNPs)
    s_plot = colSums(s.SNPs)
    m_plot = colSums(m.SNPs)
  }else{
    a = colMeans(g1.SNPs)
    b = colMeans(g9.SNPs)
    s_plot = colMeans(s.SNPs)
    m_plot = colMeans(m.SNPs)
  }
  
  if (step >1){
    x <- seq(-1*(show.window)+step/2, show.window, step)
  }else{
    x <- seq(-1*(show.window), show.window, step)
  }
  
  if(draw_plot){
    pdf(paste(name,"_AsTSS_SNPs_d=",d,"_step=",step,".pdf",sep=""), width=15, height = 10)
    par(mfcol=c(3,3))
    #par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
    #par(mgp=c(3,1,0))
    #par(cex.lab=2.2, cex.axis=2.2)  
    
    plot(x, b, col="black", xlab="distance to maxTSN", ylab="SNPs mean", main=paste(name, "all TSS d=",d," step=",step, sep=" "), type="o", ylim=c(0,max(a,b,s_plot,m_plot)), las=1)
    points(x, a,col="purple", type="o")
    legend("topleft", legend=c(paste("fdr<=0.1, n=", dim(g1)[1]), paste("fdr>0.9, n=", dim(g9)[1])),
           col=c("purple", "black"), bty = "n", lty=1, pch=1)
    
    plot(x,a-b, type="o", xlab="distance to maxTSN",ylab="substract", las=1, main="fdr0.1 - fdr0.9")#, ylim=c(-0.05,0.05))
    abline(h=0, col="green")
    
    plot(x,s_plot-m_plot, type="o", xlab="distance to maxTSN",ylab="single - multiple", las=1, main="single - multiple")#, ylim=c(-0.05,0.05))
    abline(h=0, col="green")
    abline(v=-5, col="orange")
    abline(v=5, col="orange")
    abline(v=35, col="orange")
    
    plot(x, m_plot, col="blue", xlab="distance to maxTSN", ylab="SNPs mean", main="FDR<=0.1 TSS", type="o", ylim=c(0,max(a,b,s_plot,m_plot)), las=1)
    points(x, s_plot,col="red", type="o")
    legend("topleft", legend=c(paste("Single base, n=", dim(s)[1],sep=""), paste("Multiple base, n=", dim(m)[1],sep="")),
           col=c("red", "blue"), bty = "n", lty=1, pch=1)
    

    
    plot(x, b,col="black",  xlab="distance to maxTSN", ylab="SNPs mean", main="TSS", type="o", ylim=c(0,max(a,b,s_plot,m_plot)), las=1)
    points(x, s_plot, col="red", type="o")
    
    legend("topleft", legend=c(paste("Single base, n=", dim(s)[1],sep=""), paste("fdr>0.9, n=", dim(g9)[1])),
           col=c("red", "black"), bty = "n", lty=1, pch=1)
    
    plot(x, m_plot, col="blue", xlab="distance to maxTSN", ylab="SNPs mean", main="TSS", type="o", ylim=c(0,max(a,b,s_plot,m_plot)), las=1)
    points(x, b,col="black", type="o")
    
    legend("topleft", legend=c(paste("Multiple base, n=", dim(m)[1],sep=""), paste("fdr>0.9, n=", dim(g9)[1])),
           col=c("blue", "black"), bty = "n", lty=1, pch=1)
    
    plot(x, s_plot - b,col="red", pch=19,  xlab="distance to maxTSN", ylab="SNPs mean", main="fdr 0.1 - fdr0.9 TSS", type="o", las=1)
    points(x, m_plot - b,col="blue", pch=19, type="o")
    abline(h=0, col="black") 
    #abline(v=0, col="black") 
    legend("topleft", legend=c(paste("Single base, n=", dim(s)[1],sep=""), paste("Multiple base, n=", dim(m)[1],sep="")),
           col=c("red", "blue"), bty = "n", lty=1, pch=19)
    
    plot(x, s_plot - b,col="red", pch=19, xlab="distance to maxTSN", ylab="SNPs mean", main="fdr 0.1 Single - fdr0.9 TSS", type="o", las=1)
    abline(h=0, col="green")
    
    plot(x, m_plot - b,col="blue", pch=19,  xlab="distance to maxTSN", ylab="SNPs mean", main="fdr 0.1 multiple - fdr0.9 TSS", type="o", las=1)
    abline(h=0, col="green")  
    dev.off()
    
    
  }
  return (list(x,s_plot, m_plot, dim(s), dim(m)))
  
}

Organ=c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")
s_count <- NULL
m_count <- NULL
for(t in Organ){
  #SNPsAbundanceAroundMaxTSNInTSS(maxTSN, d=50, step=1,times=1, use.log=FALSE, use.sum=FALSE, name="BN", OnlyAsTSS= FALSE)
  temp = SNPsAbundanceAroundMaxTSNInTSS(maxTSN, d=202, step=4,times=1, use.log=FALSE, use.sum=FALSE, name=t, OnlyAsTSS= TRUE)
  s_count = c(s_count, temp[[4]][1])
  m_count = c(m_count, temp [[5]][1])

  
  #SNPsAbundanceAroundMaxTSNInTSS(maxTSN, d=202, step=6,times=1, use.log=FALSE, use.sum=FALSE, name=t, OnlyAsTSS= TRUE)
  #SNPsAbundanceAroundMaxTSNInTSS(maxTSN, d=205, step=10,times=1, use.log=FALSE, use.sum=FALSE, name=t, OnlyAsTSS= TRUE)
  
  #SNPsAbundanceAroundMaxTSNInTSS(maxTSN, d=50, step=1,times=1, use.log=FALSE, use.sum=FALSE, name="BN", OnlyAsTSS= FALSE, OnlynonAsTSSwithAsMaxTSN= TRUE, p_value_cut=0.01) 
  }
df = data.frame(Organ, s_count, m_count)


# test the SNPs abundance 
# emxamine if tyhere are enrichment of SNPs in multiple base driven as.TSS in from 5 to 35bp, compare to -5 to 5 bp aournd maxTSN
p.value <- NULL
odds.ratio <- NULL
testors <- NULL
range=c(-5,5,5,35)
Organ=c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")
for(t in Organ){
  count = SNPsAbundanceAroundMaxTSNInTSS(maxTSN, d=205, step=10,times=1, use.log=FALSE, use.sum=TRUE, name=t, OnlyAsTSS= TRUE, draw_plot = FALSE)
  x = count[[1]]
  s_plot = count [[2]]
  m_plot = count [[3]]
  
  inside = which(x > range[1] & x < range[2] )
  outside = which(x > range[3] & x < range[4])
  
  testors= rbind(testors, c(sum(s_plot[inside]),sum(s_plot[outside]), sum(m_plot[inside]),sum(m_plot[outside]))) ; 
  
  testor = rbind(c(sum(s_plot[inside]),sum(s_plot[outside])),
                 + c(sum(m_plot[inside]),sum(m_plot[outside])) ); 
  cat(testor, "\n")
  f = fisher.test(testor); f
  p.value= c(p.value, f$p.value)
  odds.ratio = c(odds.ratio, f$estimate)
}
p.value
odds.ratio
adjust.p = p.adjust(p.value, method="bonferroni")
testors = data.frame(testors)
colnames(testors) = c("s1","s2","m1","m2")
data.frame(Organ, testors, odds.ratio,p.value, adjust.p )


# examine the number of as.maxTSNs in single or multiple base driven asTSS
p.value <- NULL
odds.ratio <- NULL
testors <- NULL
Organ=c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")
for(t in Organ){
  name=t
  OnlyAsTSS=TRUE
  p_value_cut=NULL
  AsTSS_maxTSN=read.table(file = paste(name,"_allReads_TSS_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue_As.maxTSNs.bed",sep=""), header = F)
  colnames(AsTSS_maxTSN)[7:10]=c("masked_p_value","masked_fdr", "unmasked_p_value","unmasked_fdr")
  colnames(AsTSS_maxTSN)[11:15]=c("chrm","chrmStart", "chrmEnd", "BinoP", "BinoP_Rfdr")
  
  if (OnlyAsTSS){  # only use as.TSS. Remove TSS with as.maxTSN
    AsTSS_maxTSN=AsTSS_maxTSN[AsTSS_maxTSN$unmasked_fdr <= 0.1 ,]
  } else if (OnlynonAsTSSwithAsMaxTSN){  # only nonAS.TSS with as.maxTSN
    AsTSS_maxTSN=AsTSS_maxTSN[AsTSS_maxTSN$unmasked_fdr > 0.1 ,]
  }
  if (is.null(p_value_cut)){
    p_value_cut = max(AsTSS_maxTSN$unmasked_p_value[AsTSS_maxTSN$unmasked_fdr<=0.1])
  }
  
  # asTSS driven by multiple base
  m = AsTSS_maxTSN[AsTSS_maxTSN$masked_p_value <= p_value_cut,c(11:16,7)]
  # asTSS driven by single base
  s = AsTSS_maxTSN[AsTSS_maxTSN$masked_p_value > p_value_cut,c(11:16,7)]
  
  # how many single base driven asTSS contains as.maxTSN?
  # how many multiple base driven asTSS contains as.maxTSN?
  cat(name, "\n")
  cat((sum(s$chrmStart != -1 & s$BinoP_Rfdr <= 0.1) / dim(s)[1])/(sum(m$chrmStart != -1 & m$BinoP_Rfdr <= 0.1) / dim(m)[1]), "\n")
  
  testor = rbind(c(sum(s$chrmStart != -1 & s$BinoP_Rfdr <= 0.1), dim(s)[1]),
                 + c(sum(m$chrmStart != -1 & m$BinoP_Rfdr <= 0.1), dim(m)[1])); 
  testors=rbind(testors, as.vector(t(testor)))
  cat(testor, "\n")
  f = fisher.test(testor); f
  p.value= c(p.value, f$p.value)
  odds.ratio = c(odds.ratio, f$estimate)
}
p.value
odds.ratio
adjust.p = p.adjust(p.value, method="bonferroni")
testors = data.frame(testors)
colnames(testors) = c("As.maxTSN in s","s","As.maxTSN in m","m")
data.frame(Organ, testors, odds.ratio,p.value, adjust.p )

if(0) {
  hist(s$BinoP_Rfdr, col="red", density = 25, breaks = seq(-1,1,0.01), right=F)
  hist(m$BinoP_Rfdr, add=T, col="blue", breaks = seq(-1,1,0.01), right=F)
  hist(s$BinoP_Rfdr, add=T, col="red", density = 25, breaks = seq(-1,1,0.01), right=F)
  legend("topright", 
         legend = c("multiple", "single"),
         #pch=c(15,15),
         cex=2, 
         lty=c(0,0),
         #bty="n",
         lwd=1.5, 
         density=c(10000,25),
         angle=c(180,45),
         #angle=45,
         fill=c("blue","red")
         , bty = "n"
  )
}


# seperate as.maxTSN in three group
# inside s, single base drive asTSS
# inside m, multiple base driven asTSS
# others
SNPsAbundanceAroundAsMaxTSNInTSS <-function(d=50, step=1,times=1, use.log=FALSE, use.sum=FALSE, name="", OnlyAsTSS= TRUE, OnlynonAsTSSwithAsMaxTSN= FALSE, p_value_cut=NULL , draw_plot=TRUE){ 

#d=202
  #step=4
  # only contain maxTSN with SNPs within 20bp downstream
  maxTSN_20bpSNP=read.table(file = paste(name,"_allReads_TSS_5mat5pat_uniq_WithAsMaxTSNunionAsTSSfdr0.1_maskedVSunmasked_pValue_As.maxTSNs_ba.bed",sep=""), header = F)
  colnames(maxTSN_20bpSNP)[13:16]=c("masked_p_value","masked_fdr", "unmasked_p_value","unmasked_fdr")
  colnames(maxTSN_20bpSNP)[1:5]=c("chrm","chrmStart", "chrmEnd","p_value", "Rfdr") 
  
  NonAsmaxTSN = maxTSN_20bpSNP[maxTSN_20bpSNP$Rfdr > 0.9,]
  AsmaxTSN = maxTSN_20bpSNP[maxTSN_20bpSNP$Rfdr <=0.1,]

  withoutTSSAStest = AsmaxTSN[AsmaxTSN$V17==0,1:6]
  withTSSAStest = AsmaxTSN[AsmaxTSN$V17==1,]
  
  withTSSAStest$masked_p_value=as.numeric(as.character(withTSSAStest$masked_p_value))
  withTSSAStest$unmasked_p_value=as.numeric(as.character(withTSSAStest$unmasked_p_value))
  withTSSAStest$unmasked_fdr=as.numeric(as.character(withTSSAStest$unmasked_fdr))
  
  nonAsTSS = withTSSAStest[withTSSAStest$unmasked_fdr>0.1,1:6]
  inAsTSS = withTSSAStest[withTSSAStest$unmasked_fdr<=0.1,]
  
  p_value_cut = max(inAsTSS$unmasked_p_value)
  # 
  # asTSS driven by multiple base
  m = inAsTSS [(inAsTSS $masked_p_value <= p_value_cut),]
  # asTSS driven by single base
  s = inAsTSS[(inAsTSS$masked_p_value > p_value_cut),]
  

  
  cat(name, "\n")
  cat ("dim(NonAsmaxTSN)[1]:" )
  cat (dim(NonAsmaxTSN)[1],  "\n")
  cat ("dim(AsmaxTSN)[1]:" )
  cat (dim(AsmaxTSN)[1],  "\n")
  cat ("dim(withoutTSSAStest)[1]:" )
  cat (dim(withoutTSSAStest)[1],  "\n")
  cat ("dim(nonAsTSS)[1]:" )
  cat (dim(nonAsTSS)[1],  "\n")
  cat ("dim(inAsTSS)[1]:" )
  cat (dim(inAsTSS)[1],  "\n")
  
  cat ("dim(m)[1], max(m$masked_p_value), min(m$masked_p_value)", "\n")
  cat (dim(m)[1], max(m$masked_p_value), min(m$masked_p_value), "\n")
  cat("  dim(s)[1]; max(s$masked_p_value); min(s$masked_p_value)", "\n")
  cat(dim(s)[1], max(s$masked_p_value), min(s$masked_p_value), "\n")
  
  NonAsmaxTSN$chrmStart = NonAsmaxTSN$chrmStart-d
  NonAsmaxTSN$chrmEnd = NonAsmaxTSN$chrmEnd+d
  AsmaxTSN$chrmStart = AsmaxTSN$chrmStart-d
  AsmaxTSN$chrmEnd = AsmaxTSN$chrmEnd+d
  withoutTSSAStest$chrmStart = withoutTSSAStest$chrmStart-d
  withoutTSSAStest$chrmEnd = withoutTSSAStest$chrmEnd+d
  nonAsTSS$chrmStart = nonAsTSS$chrmStart-d
  nonAsTSS$chrmEnd = nonAsTSS$chrmEnd+d
  inAsTSS$chrmStart = inAsTSS$chrmStart-d
  inAsTSS$chrmEnd = inAsTSS$chrmEnd+d
  s$chrmStart = s$chrmStart-d
  s$chrmEnd = s$chrmEnd+d
  m$chrmStart = m$chrmStart-d
  m$chrmEnd = m$chrmEnd+d
  

  
  # get the sum of SNPs in each bin (size = step)
  NonAsmaxTSN.SNPs <- read_read_mat_SNPs (SNP.bw, NonAsmaxTSN[,c(1:6)], step, times=times, use.log=use.log)
  AsmaxTSN.SNPs <- read_read_mat_SNPs (SNP.bw, AsmaxTSN[,c(1:6)], step, times=times, use.log=use.log)
  withoutTSSAStest.SNPs <- read_read_mat_SNPs (SNP.bw, withoutTSSAStest[,c(1:6)], step, times=times, use.log=use.log)
  nonAsTSS.SNPs <- read_read_mat_SNPs (SNP.bw, nonAsTSS[,c(1:6)], step, times=times, use.log=use.log)
  inAsTSS.SNPs <- read_read_mat_SNPs (SNP.bw, inAsTSS[,c(1:6)], step, times=times, use.log=use.log)
  
  
  s.SNPs <- read_read_mat_SNPs (SNP.bw, s[,c(1:6)], step, times=times, use.log=use.log)
  m.SNPs <- read_read_mat_SNPs (SNP.bw, m[,c(1:6)], step, times=times, use.log=use.log)
 
  show.window <- d
  if (use.sum){
    NonAsmaxTSN_plot = colSums(NonAsmaxTSN.SNPs)
    AsmaxTSN_plot = colSums(AsmaxTSN.SNPs)
    withoutTSSAStest_plot = colSums(withoutTSSAStest.SNPs)
    nonAsTSS_plot = colSums(nonAsTSS.SNPs)
    inAsTSS_plot = colSums(inAsTSS.SNPs)
    s_plot = colSums(s.SNPs)
    m_plot = colSums(m.SNPs)
  }else{
    NonAsmaxTSN_plot = colMeans(NonAsmaxTSN.SNPs)
    AsmaxTSN_plot = colMeans(AsmaxTSN.SNPs)
    withoutTSSAStest_plot = colMeans(withoutTSSAStest.SNPs)
    nonAsTSS_plot = colMeans(nonAsTSS.SNPs)
    inAsTSS_plot = colMeans(inAsTSS.SNPs)
    s_plot = colMeans(s.SNPs)
    m_plot = colMeans(m.SNPs)
  }
  
  if (step >1){
    x <- seq(-1*(show.window)+step/2, show.window, step)
  }else{
    x <- seq(-1*(show.window), show.window, step)
  }
  
  if(draw_plot){
    pdf(paste(name,"_AsMaxTSN_TSS_SNPs_d=",d,"_step=",step,"-1.pdf",sep=""), width=10, height = 10)
    #par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
    #par(mgp=c(3,1,0))
    #par(cex.lab=2.2, cex.axis=2.2)  

    plot(x, AsmaxTSN_plot,col="green", xlab="distance to ASmaxTSN", ylab="SNPs mean", main=paste(name, "all TSS d=",d," step=",step, sep=" "), type="o", ylim=c(0,max(AsmaxTSN_plot,inAsTSS_plot,s_plot,m_plot)), las=1)
    points(x, nonAsTSS_plot,col="black", type="o", pch=19)
    points(x, inAsTSS_plot,col="purple", type="o", pch=19) #,bg="black")
    points(x, withoutTSSAStest_plot,col="brown", type="o")
    points(x, AsmaxTSN_plot,col="green", type="o")
    
    legend("topleft", legend=c(paste("As.maxTSN All, n=", dim(AsmaxTSN)[1]), 
                                paste("As.maxTSN in asTSS, n=", dim(inAsTSS)[1]), 
                                paste("As.maxTSN in non.asTSS, n=", dim(nonAsTSS)[1]), 
                                paste("As.maxTSN in TSS without AS test, n=", dim(withoutTSSAStest)[1])),
           col=c("green", "purple","black", "brown"), bty = "n", lty=1, pch=c(1,19,19,1))
    dev.off()
    
    
    pdf(paste(name,"_AsMaxTSN_TSS_SNPs_d=",d,"_step=",step,".pdf",sep=""), width=10, height = 10)
    par(mfrow=c(4,2))
    
    plot(x, AsmaxTSN_plot,col="green", xlab="distance to maxTSN", ylab="SNPs mean", main=paste(name, "all maxTSN_20bp d=",d," step=",step, sep=" "), type="o", ylim=c(0,max(AsmaxTSN_plot,inAsTSS_plot,s_plot,m_plot)), las=1)
    points(x, NonAsmaxTSN_plot,col="black", type="o", pch=17, bg="black")    
    legend("topleft", legend=c(paste("As.maxTSN (fdr<=0.1), n=", dim(AsmaxTSN)[1]), 
                               paste("Non.As.maxTSN (fdr>0.9), n=", dim(NonAsmaxTSN)[1])), 
           col=c("green", "black"), bty = "n", lty=1, pch=c(1,17))
    
    plot(x,AsmaxTSN_plot - NonAsmaxTSN_plot, type="o", xlab="distance to maxTSN",ylab="substract", las=1, main="asMaxTSN(fdr0.1) - non.as.MaxTSN(fdr0.9)")#, ylim=c(-0.05,0.05))
    abline(h=0, col="green")
    
    plot(x, AsmaxTSN_plot,col="green", xlab="distance to ASmaxTSN", ylab="SNPs mean", main=paste(name, "all AsMaxTSN d=",d," step=",step, sep=" "), type="o", ylim=c(0,max(AsmaxTSN_plot,inAsTSS_plot,s_plot,m_plot)), las=1)
    points(x, nonAsTSS_plot,col="black", type="o", pch=19)
    points(x, inAsTSS_plot,col="purple", type="o", pch=19) #,bg="black")
    points(x, withoutTSSAStest_plot,col="brown", type="o")
    points(x, AsmaxTSN_plot,col="green", type="o")
    
    legend("topleft", legend=c(paste("As.maxTSN All, n=", dim(AsmaxTSN)[1]), 
                               paste("As.maxTSN in asTSS, n=", dim(inAsTSS)[1]), 
                               paste("As.maxTSN in non.asTSS, n=", dim(nonAsTSS)[1]), 
                               paste("As.maxTSN in TSS without AS test, n=", dim(withoutTSSAStest)[1])),
           col=c("green", "purple","black", "brown"), bty = "n", lty=1, pch=c(1,19,19,1))
    
    plot(x,inAsTSS_plot- nonAsTSS_plot, type="o", xlab="distance to ASmaxTSN",ylab="substract", las=1, main="asMaxTSN in asTSS(fdr0.1) - asMaxTSN in nonasTSS(fdr0.9)")#, ylim=c(-0.05,0.05))
    abline(h=0, col="green")
    
    plot(x, m_plot,col="blue", xlab="distance to ASmaxTSN", ylab="SNPs mean", main="asMaxTSN in asTSS(fdr<=0.1)", type="o", ylim=c(0,max(AsmaxTSN_plot,inAsTSS_plot,s_plot,m_plot)), las=1)
    points(x, s_plot,col="red", type="o")
    legend("topleft", legend=c(paste("AsMaxTSN in Single base, n=", dim(s)[1],sep=""), paste("AsMaxTSN in Multiple base, n=", dim(m)[1],sep="")),
           col=c("red", "blue"), bty = "n", lty=1, pch=1)
    
    plot(x,s_plot-m_plot, type="o", xlab="distance to ASmaxTSN",ylab="single - multiple", las=1, main="single - multiple")#, ylim=c(-0.05,0.05))
    abline(h=0, col="green")
    abline(v=-5, col="orange")
    abline(v=5, col="orange")
    abline(v=35, col="orange")
    
    if(0){
    plot(x, inAsTSS_plot, pch=19, col="purple",  xlab="distance to ASmaxTSN", ylab="SNPs mean", main="TSS", type="o", ylim=c(0,max(AsmaxTSN_plot,inAsTSS_plot,s_plot,m_plot)), las=1)
    points(x, s_plot, col="red", type="o")
    points(x, m_plot, col="blue", type="o")
    
    legend("topleft", legend=c(paste("AsTSS Single base, n=", dim(s)[1],sep=""), paste("AsTSS Multiple base, n=", dim(m)[1],sep=""), paste("All AsTSS fdr<0.1, n=", dim(inAsTSS)[1])),
           col=c("red", "blue", "purple"), bty = "n", lty=1, pch=c(1,1,19))
    }

    plot(x, nonAsTSS_plot,col="black",  pch=19,
         xlab="distance to ASmaxTSN", ylab="SNPs mean", main="AsMaxTSN in TSS", type="o", ylim=c(0,max(AsmaxTSN_plot,inAsTSS_plot,s_plot,m_plot)), las=1)
    points(x, s_plot, col="red", type="o")
    legend("topleft", legend=c(paste("AsMaxTSN in SingleD AsTSS , n=", dim(s)[1],sep=""), paste("AsMaxTSN in nonAs TSS fdr>0.9, n=", dim(nonAsTSS)[1])),
           col=c("red", "black"), bty = "n", lty=1, pch=c(1,19))
    
   
    plot(x, nonAsTSS_plot,col="black",  pch=19,
         xlab="distance to ASmaxTSN", ylab="SNPs mean", main="TSS", type="o", ylim=c(0,max(AsmaxTSN_plot,inAsTSS_plot,s_plot,m_plot)), las=1)
    points(x, m_plot, col="blue", type="o")
    legend("topleft", legend=c(paste("AsMaxTSN in MultipleD AsTSS, n=", dim(m)[1],sep=""), paste("AsMaxTSN in nonAs TSS fdr>0.9, n=", dim(nonAsTSS)[1])),
           col=c("blue", "black"), bty = "n", lty=1, pch=c(1,19))
    
    dev.off()
  }

  
}

for(t in c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")){
#  t="BN"
  SNPsAbundanceAroundAsMaxTSNInTSS (d=50, step=1,times=1, use.log=FALSE, use.sum=FALSE, name=t, OnlyAsTSS= TRUE, OnlynonAsTSSwithAsMaxTSN= FALSE, p_value_cut=NULL , draw_plot=TRUE)
  SNPsAbundanceAroundAsMaxTSNInTSS (d=51, step=2,times=1, use.log=FALSE, use.sum=FALSE, name=t, OnlyAsTSS= TRUE, OnlynonAsTSSwithAsMaxTSN= FALSE, p_value_cut=NULL , draw_plot=TRUE)
  SNPsAbundanceAroundAsMaxTSNInTSS (d=52, step=4,times=1, use.log=FALSE, use.sum=FALSE, name=t, OnlyAsTSS= TRUE, OnlynonAsTSSwithAsMaxTSN= FALSE, p_value_cut=NULL , draw_plot=TRUE)
  SNPsAbundanceAroundAsMaxTSNInTSS (d=52, step=6,times=1, use.log=FALSE, use.sum=FALSE, name=t, OnlyAsTSS= TRUE, OnlynonAsTSSwithAsMaxTSN= FALSE, p_value_cut=NULL , draw_plot=TRUE)
  
}

