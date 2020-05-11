# Rscript getSNPsAbundance.R

require(bigWig)


### SNP locations
#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/Initiation")
#source("/Users/shaopei/Box\ Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/heatmap/heatmaps.R")

SNP.bw <- "/Volumes/SPC_SD/KD_IGV/P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"
#SNP.bw <- "P.CAST_M.B6_indelsNsnps_CAST.bam.snp.unfiltered_plus.bw"

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
  
  plot(x, b, col="blue", xlab="distance to maxSNP", ylab="SNPs mean", main=name, type="o", ylim=c(0,max(a,b)), las=1)
  points(x, a,col="red", type="o")
  legend("topright", legend=c("fdr<=0.1", "fdr>0.9"),
         col=c("red", "blue"), bty = "n", lty=1, pch=1)
  plot(x,a-b, type="o", xlab="distance to maxSNP",ylab="substract", las=1, main="substract")#, ylim=c(0,max(a,b)))
  dev.off()

}

for(t in c("BN","LV","HT", "SK", "KD", "SP", "GI", "ST")){
df=read.table(file = paste(t,"_allReads_TSS_maxTSNs_SNPs20bp_binomtest_Rfdr1.bed",sep = ""), header = T)
SNPsAbundanceAroundMaxTSN(df, name=t)
}
