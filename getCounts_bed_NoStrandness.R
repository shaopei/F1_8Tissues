##
## getCounts_bed.R - Counts reads in each region of the provided bed file using bigWig file as input
args=(commandArgs(TRUE))
setwd(args[1])
df_fp=args[2] #"T8_2Strand_p0.05_effect_strain.bed_cluster_organSpecific"

#setwd("~/Box Sync/Danko_lab_work/F1_8Tissues")
#df_fp = "T8_2Strand_p0.05_effect_strain.bed_cluster_organSpecific"

require(bigWig)

bodies <- read.table(df_fp, header=F)

countBigWig <- function(prefix, bed, rpkm=rpkm, path="./") {
 pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))

 counts_pl <- bed.region.bpQuery.bigWig(pl, bed, abs.value = TRUE)
 counts_mn <- bed.region.bpQuery.bigWig(mn, bed, abs.value = TRUE)
 counts = counts_pl + counts_mn
 
        if(rpkm==TRUE) {
                counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }

 return(counts)
}

stage     <- c("BN", "HT", "SK", "SP", "KD", "LV", "GI" ,"ST")
replicate <- c("all")

filenames=c()
for (s in stage){
  filenames <- c(filenames, paste(s, replicate, sep="_"))
}
## Gets counts
counts <- NULL
for(f in filenames) {
  counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE))
}
colnames(counts) <- filenames
counts <- cbind(counts, bodies)
save.image("data-counts.RData")
#remove(counts)


## Gets RPKMs
rpkm <- NULL
for(f in filenames) {
    rpkm <- cbind(rpkm, countBigWig(f, bodies, rpkm=TRUE))
}
colnames(rpkm) <- stage
rpkm <- cbind(rpkm, bodies)

rpkm$B=0 # rpkm of the organ specific domain in Biased organ)
rpkm$H=0 # rpkm of the non-biased organ with highest expression level
for (i in 1:dim(rpkm)[1]){
  rpkm$B[i]=rpkm[i,grep(rpkm$V4[i],colnames(rpkm))]
  rpkm$H[i]=max(rpkm[i,grep(rpkm$V4[i],colnames(rpkm)[1:8], invert = T )])
  rpkm$NonBH[i]=colnames(rpkm)[grep(rpkm$H[i],rpkm[i,1:8])]
}

save.image("data-rpkms.RData")
load("data-rpkms.RData")
rpkm$NonBH_B_ratio = rpkm$H/rpkm$B
pdf("allellicBias_expression_level1.pdf")
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)

hist(log2(rpkm$NonBH_B_ratio), xlab="log2(rpkm non-biased organ with highest expression level/Biased organ)")
dev.off()


