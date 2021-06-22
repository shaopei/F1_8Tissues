## cd /workdir/sc2457/F1_Tissues/map2ref_1bpbed_map5_MultiBaseRunOn/map2ref_1bpbed_map5
# R
## getCounts.R - Counts reads in each gene.
require(bigWig)

countBigWig <- function(prefix, bed, rpkm=FALSE, path="./") {
 pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))

 counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
        if(rpkm==TRUE) {
                counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }

 return(counts)
}

bodies <- read.table("/workdir/sc2457/F1_Tissues/RNA-seq/STAR_BN/gencode.vM25.annotation.gene.bed", header=F)
bodies <- bodies[bodies$V1 != "chrX",]
bodies <- bodies[bodies$V1 != "chrY",]
bodies <- bodies[bodies$V1 != "chrM",]

filenames <- c("LV_map2ref_1bpbed_map5", "LV_map2ref_1bpbed_map5_B6", "LV_map2ref_1bpbed_map5_CAST")

## Gets counts
counts <- NULL
for(f in filenames) {
	counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE))
    #pause_counts <- cbind(pause_counts, countBigWig(f, pause, rpkm=FALSE))
	#postcps_counts <- cbind(postcps_counts, countBigWig(f, postcps, rpkm=FALSE)) 
}

colnames(counts) <- c("All.proseq", "B6.proseq", "CAST.proseq")
counts <- cbind(bodies, counts)
colnames(counts)[4]="geneID"
colnames(counts)[5]="geneName"

#exon.filenames <- c("LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.exon.read.count",
exon.filenames <- c(	"LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.nsorted.exon.read.count", 
	"LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast.nsorted.exon.read.count")
exon.counts <- NULL
for(f in exon.filenames) {
	r <- read.table(paste("/workdir/sc2457/F1_Tissues/RNA-seq/STAR_LV/",f,sep=""), header=F)
	#merge.counts <- merge(rna.counts, r, by="V1", all=TRUE) 
	exon.counts <- cbind(exon.counts, r[,2])
}
exon.counts <- cbind.data.frame(r$V1, exon.counts)
#colnames(exon.counts) <- c("geneID","All.exon", "B6.exon", "CAST.exon")
colnames(exon.counts) <- c("geneID","B6.exon", "CAST.exon")

if(0){
rna.filenames <- c("LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.gene.read.count",
	"LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.nsorted.gene.read.count", 
	"LV_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast.nsorted.gene.read.count")
rna.counts <- NULL
for(f in rna.filenames) {
	r <- read.table(paste("/workdir/sc2457/F1_Tissues/RNA-seq/STAR_LV/",f,sep=""), header=F)
	#merge.counts <- merge(rna.counts, r, by="V1", all=TRUE) 
	rna.counts <- cbind(rna.counts, r[,2])
}
rna.counts <- cbind.data.frame(r$V1, rna.counts)
colnames(rna.counts) <- c("geneID","All.rna", "B6.rna", "CAST.rna")
}

#merge.counts <- merge(counts, rna.counts, by="geneID", all.x=TRUE) 
merge.counts <- merge(counts, exon.counts, by="geneID", all.x=TRUE) 

save.image("LV_gencode.vM25.annotation.gene-Readcounts.RData")

setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/RNA-seq")
load("LV_gencode.vM25.annotation.gene-Readcounts.RData")

counts = merge.counts
dim(counts)
counts = counts[counts$B6.proseq >= 10,]
counts = counts[counts$CAST.proseq >= 10,]
counts$sB6 <- counts$B6.exon/counts$B6.proseq
counts$sCAST <- counts$CAST.exon/counts$CAST.proseq
dim(counts)

plot(counts$sB6, counts$sCAST,
     xlim=c(0,200), ylim = c(0,200))


geneID_withATwindow <- read.table("LV_geneID_withATwindow", header=F)
geneID_withATwindow$ATwindow = TRUE
geneID_withATwindow_withRNAAlleleHMMBlocks <- read.table("LV_geneID_withATwindow.with.nearby.RNA.AlleleHMM.blocks", header=F)
geneID_withATwindow_withRNAAlleleHMMBlocks$AlleleHMM <- 1

target = merge(geneID_withATwindow, geneID_withATwindow_withRNAAlleleHMMBlocks, by = "V1", all.x = T)
dim(target)
colnames(target)[1]="geneID"
target$AlleleHMM[is.na(target$AlleleHMM)] = 0
target = merge(target, counts, by = "geneID", all.x = T)
target$deltaS = abs(target$sB6-target$sCAST)
dim(counts)
dim(target)
target=target[!(is.na(target$V1)),]
dim(target)

LV_t=target
write.table(LV_t, file="LV_geneWithInATwindows_stability.txt", quote = F, sep="\t", col.names = T, row.names = F)
#hist(target$sB6, breaks = seq(0,14,0.01),col="blue")
#hist(target$sCAST, breaks = seq(0,14,0.01), col="red")

sum(target$AlleleHMM==0)
sum(target$AlleleHMM==1)
pdf("LV_Allelic_Termination_Gene_Stablility.pdf", width=7, height = 3.5, useDingbats=FALSE)
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
plot(ecdf(target$deltaS[target$AlleleHMM==0]),
     #pch=10, 
     col="dark blue",
     xlim=c(0,2),
     ylim=c(0.9,1),
     xlab="|Stability B6 - Stability CAST|",
     ylab="CDF",
     las=1,cex=1,
     main=""
)
lines(ecdf(target$deltaS[target$AlleleHMM==1]), cex=1, col="dark orange")

legend("right", 
       legend = c( paste("Genes without allelic difference in mature mRNA, n = ", sum(target$AlleleHMM==0), sep=""), 
                   #paste("Without Indel, n=", dim(sub_df_no_indel)[1], sep=""),
                   paste("Genes with allelic difference in mature mRNA, n=", sum(target$AlleleHMM==1), sep="")),
       title = ,
       #pch=c(15,15),
       #lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       #fill=c("blue","dark organe","dark green")
       col=c("dark blue","dark orange"),
       pch=c(19,19),
       bty = "n"
)
dev.off()

# KS test 
#ks.test(target$deltaS[target$AlleleHMM==0] ,target$deltaS[target$AlleleHMM==1], alternative = "two.sided")
ks.test(target$deltaS[target$AlleleHMM==0] ,target$deltaS[target$AlleleHMM==1], alternative = "greater")
ks.test(target$deltaS[target$AlleleHMM==0] ,target$deltaS[target$AlleleHMM==1], alternative = "less")
