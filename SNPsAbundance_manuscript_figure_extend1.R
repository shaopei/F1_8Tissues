d=50; step=1;times=1; use.log=FALSE; use.sum=FALSE; name=""; 
OnlyAsTSS= FALSE; OnlynonAsTSSwithAsMaxTSN= FALSE; p_value_cut=NULL ;
draw_plot=TRUE;
s_l=NULL; m_l=NULL


name="BN"
d=50; OnlyAsTSS= TRUE


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

#hist(rowSums(s.SNPs), breaks = seq(0,20,1))
#hist(rowSums(m.SNPs), breaks = seq(0,20,1))
mean(rowSums(s.SNPs))
mean(rowSums(m.SNPs))
ks.test(rowSums(s.SNPs), rowSums(m.SNPs), alternative = "less")
plot(ecdf(rowSums(m.SNPs)), col="blue", 
     xlab= paste("SNPs count per site (d= ", d, " bp)", sep=""),
     ylab="density",
     las=1,
     main=name
     )
lines(ecdf(rowSums(s.SNPs)), col="red")
hist(rowSums(s.SNPs), breaks = seq(0,20,1), freq = F, add=T, col="red", density = 45, angle = -45)
hist(rowSums(m.SNPs), breaks = seq(0,20,1), freq = F, add=T, col = "blue", density = 45)

legend("right", 
       legend = c( paste("single SNPs per row, n=", dim(s)[1], sep=""),
         paste("multiple SNPs per row, n=", dim(m)[1], sep="")),
       title = ,
       #pch=c(15,15),
       #cex=2, 
       lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       density=45,
       angle=c(-45, 45),
       #angle=45,
       fill=c("red","blue")
       , bty = "n"
)


SNPsAbundanceAroundMaxTSNInTSS <-function(d=50, step=1,times=1, use.log=FALSE, use.sum=FALSE, name="", 
                                          OnlyAsTSS= FALSE, OnlynonAsTSSwithAsMaxTSN= FALSE, p_value_cut=NULL ,
                                          draw_plot=TRUE,
                                          s_l=NULL, m_l=NULL){ 
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
  
    hist(rowSums(s.SNPs), breaks = seq(0,20,1))
    hist(rowSums(m.SNPs), breaks = seq(0,20,1))
    ks.test(rowSums(s.SNPs), rowSums(m.SNPs), alternative = "greater")
  
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
    if(0){
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
    }
    pdf(paste(name,"_AsTSS_SNPs_d=",d,"_step=",step,"-S.pdf",sep=""), width=9.23, height = 6.52)
    par(mfcol=c(2,1))
    par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
    par(mgp=c(3,1,0))
    par(cex.lab=2.2, cex.axis=2.2)
    #par(cex=1.5)
    pch_u=15
    plot(x, b,col="black",  xlab="distance to maxTSN", ylab="SNPs mean", 
         main=name, type="o", ylim=c(0,max(a,b,s_plot,m_plot)), pch=pch_u, las=1, frame=FALSE)
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
    if(0){
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
    
  }
  return (list(x,s_plot, m_plot, dim(s), dim(m), b, dim(g9)))
  
}
