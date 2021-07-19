##
## getCounts.R - Counts reads in each gene.
file_dir="~/Box Sync/BN_IGV/"
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/")

organ="BN"
require(bigWig)

tus_g9=read.table(file = paste(organ,"_initiation_g9_NoMS_tunits.bed",sep=""), header=FALSE)
tus_g9$V12 = "g9"
tus_m=read.table(file = paste(organ,"_initiation_m_tunits.bed",sep=""), header=FALSE)
tus_m$V12 = "m"
tus_s=read.table(file = paste(organ,"_initiation_s_tunits.bed",sep=""), header=FALSE)
tus_s$V12 = "s"

tus <- rbind( tus_m, tus_s, tus_g9 )
dim(tus)
#tus <- read.table(paste("./tunit_protein_coding/",organ, "_all_h5.preds.full_inProtein_coding.bed", sep=""), header=FALSE)
colnames(tus)[9] = "TXSTART"
colnames(tus)[10] = "TXEND"
colnames(tus)[12] = "group"
colnames(tus)[13] = "TXSTRAND"
tus <- tus[(tus$TXEND-tus$TXSTART)>500,]
dim(tus)
# body exclude the first 250bp of transcript
bodies <- tus[,8:13]
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+250
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-250

countBigWig <- function(prefix, bed, rpkm=FALSE, path="./") {
 pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))

 counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
        if(rpkm==TRUE) {
                counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }

 return(counts)
}

b6.bw=paste(organ, "_map2ref_1bpbed_map5_B6", sep="")
cast.bw=paste(organ, "_map2ref_1bpbed_map5_CAST", sep="")
filenames <- c(b6.bw, cast.bw)
## Gets counts
counts <- NULL
for(f in filenames) {
	counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE, path=file_dir))
}

tus$B6_counts = counts[,1]
tus$CAST_counts = counts[,2]
tus$change=log2((tus$B6_counts+1)/(tus$CAST_counts+1))
library("vioplot")
vioplot(log2((tus$B6_counts+1)/(tus$CAST_counts+1)) ~ tus$group,
        main=organ)
abline(h=0, col="red")
boxplot(log2((tus$B6_counts+1)/(tus$CAST_counts+1)) ~ tus$group,
        main=organ)
abline(h=0, col="red")

cat("g9, N=" , sum(tus$group=="g9"))
cat("m, N=" , sum(tus$group=="m"))
cat("s, N=" , sum(tus$group=="s"))

change_ratio=2
testor <-  matrix(c(sum(abs(tus$change)>change_ratio & tus$group=="m"), sum(tus$group=="m"),
                    sum(abs(tus$change)>change_ratio & tus$group=="g9"),sum(tus$group=="g9")),
                  nrow = 2,
                  dimnames = list(TSS = c("change>1", "all"),
                                  KS.Test = c("m", "g9"))); testor

testor <-  matrix(c(sum(abs(tus$change)>change_ratio & tus$group=="s"), sum(tus$group=="s"),
                    sum(abs(tus$change)>change_ratio & tus$group=="g9"),sum(tus$group=="g9")),
                  nrow = 2,
                  dimnames = list(TSS = c("change>1", "all"),
                                  KS.Test = c("s", "g9"))); testor
change_ratio=2
testor <-  matrix(c(sum(abs(tus$change)>change_ratio & (tus$group=="s"|tus$group=="m")), sum(tus$group=="s"|tus$group=="m"),
                    sum(abs(tus$change)>change_ratio & tus$group=="g9"),sum(tus$group=="g9")),
                  nrow = 2,
                  dimnames = list(log2.change = c(paste("change>",change_ratio, sep=""), "all"),
                                  KS.Test = c("s+m", "g9"))); testor

f = fisher.test(testor, alternative = "g" ); f
#p.value= c(p.value, f$p.value)
#odds.ratio = c(odds.ratio, f$estimate)


### pause ###

setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/pause/")
file_dir="~/Box Sync/KD_IGV/"
organ="SK"
require(bigWig)
#body="_dREG_5mat5pat_uniq_pValue"  # use dREG KS test
#body9 = body
#body1 = body
body9="_allelicPauseWithSharedMaxTSN_shareMaxpauses"  # use pause from sharedMaxTSN KS test
body1="_allelicPauseWithSharedMaxTSN_distinctMaxpauses"  # use pause from sharedMaxTSN KS test

tus_g9=read.table(file = paste(organ, body9, "_fdr0.9_tunits.bed",sep=""), header=FALSE)
tus_g9$V14 = "g9"
tus_g1=read.table(file = paste(organ, body1, "_fdr0.1_tunits.bed",sep=""), header=FALSE)
tus_g1$V14 = "g1"
tus <- rbind( tus_g1,tus_g9 )
dim(tus)
#tus <- read.table(paste("./tunit_protein_coding/",organ, "_all_h5.preds.full_inProtein_coding.bed", sep=""), header=FALSE)
colnames(tus)[7] = "earlyPause"
colnames(tus)[8] = "latePause"
colnames(tus)[9] = "earlyPauseAllele"
colnames(tus)[11] = "TXSTART"
colnames(tus)[12] = "TXEND"
colnames(tus)[14] = "group"
colnames(tus)[15] = "TXSTRAND"

tus <- tus[(tus$TXEND-tus$TXSTART)>500,]
dim(tus)
# body exclude the first 250bp of transcript
bodies <- tus[,10:15]
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+250
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-250

countBigWig <- function(prefix, bed, rpkm=FALSE, path="./") {
    pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
    mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))
    
    counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
    if(rpkm==TRUE) {
        counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
    }
    
    return(counts)
}

b6.bw=paste(organ, ".mat.map2ref.1bp", sep="")
cast.bw=paste(organ, ".pat.map2ref.1bp", sep="")
filenames <- c(b6.bw, cast.bw)
## Gets counts
counts <- NULL
for(f in filenames) {
    counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE, path=file_dir))
}


tus$B6_counts = counts[,1]
tus$CAST_counts = counts[,2]
tus$log2_B6_CAST_ratio=log2((tus$B6_counts+1)/(tus$CAST_counts+1))
tus$log2_early_late_ratio=log2((tus$B6_counts+1)/(tus$CAST_counts+1))
tus$log2_early_late_ratio[tus$earlyPauseAllele=="CAST"] = log2((tus$CAST_counts+1)/(tus$B6_counts+1))[tus$earlyPauseAllele=="CAST"]

library("vioplot")
vioplot(log2((tus$B6_counts+1)/(tus$CAST_counts+1)) ~ tus$group,
        main=organ)
abline(h=0, col="red")
boxplot(log2((tus$B6_counts+1)/(tus$CAST_counts+1)) ~ tus$group,
        main=organ)
abline(h=0, col="red")

vioplot(tus$log2_early_late_ratio ~ tus$group,
        main=organ)
abline(h=0, col="red")
boxplot(tus$log2_early_late_ratio ~ tus$group,
        main=organ)
abline(h=0, col="red")

cat("g9, N=" , sum(tus$group=="g9"))
cat("g1, N=" , sum(tus$group=="g1"))


change_ratio=2
testor <-  matrix(c(sum(abs(tus$change)>=change_ratio & tus$group=="g1"), sum(tus$group=="g1"),
                    sum(abs(tus$change)>=change_ratio & tus$group=="g9"),sum(tus$group=="g9")),
                  nrow = 2,
                  dimnames = list(log2.change = c(paste("change>",change_ratio, sep=""), "all"),
                                  KS.Test = c("g1", "g9"))); testor

f = fisher.test(testor, alternative = "g" ); f


# AT, Allelic termination 
file_dir="~/Box Sync/BN_IGV/"
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/AT_AllelicTermination/")
organ="BN"
#tus=read.table("BN_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed") 
tus_noAT=read.table(file = paste(organ,"_all_h5.preds.full_inProtein_coding_withoutATwindow_f0.5F0.8gencode.vM25.annotation.gene.bed",sep=""), header=FALSE)
# tus_noAT col1-6 tunits 7-12 gene col13 overlap base counts
#tus_noAT = tus_noAT[,1:6] 
tus_noAT$AT_long = "NoAT"
tus_noAT$group = "NoAT"
tus_noAT=unique(tus_noAT)
dim(tus_noAT)

tus_withAT=read.table(file = paste(organ,"_AT_4tunitIntersectNativeHMM_intersectRegion_strain_f0.5F0.8gencode.vM25.annotation.gene.bed",sep=""), header=FALSE)
tus_withAT = tus_withAT[,c(13:18,20:26,19)]  #tunits
tus_withAT$group = "AT"
colnames(tus_withAT)[1:15] = colnames(tus_noAT)[1:15]
tus_withAT=unique(tus_withAT)
dim(tus_withAT)

tus=rbind(tus_noAT, tus_withAT)
dim(tus)
#View(tus)
tus=unique(tus)  #remove duplicates
dim(tus)

colnames(tus)[2] = "TXSTART"
colnames(tus)[3] = "TXEND"
colnames(tus)[6] = "TXSTRAND"
colnames(tus)[10] = "geneID"
tus <- tus[(tus$TXEND-tus$TXSTART)>500,]

dim(tus)
# body exclude the first 250bp of transcript
bodies <- tus[,1:6]
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+250
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-250

b6.bw=paste(organ, "_map2ref_1bpbed_map5_B6", sep="")
cast.bw=paste(organ, "_map2ref_1bpbed_map5_CAST", sep="")
filenames <- c(b6.bw, cast.bw)
## Gets counts
counts <- NULL
for(f in filenames) {
    counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE, path=file_dir))
}


tus$B6_counts = counts[,1]
tus$CAST_counts = counts[,2]
tus$log2_B6_CAST_ratio=log2((tus$B6_counts+1)/(tus$CAST_counts+1))
tus$log2_long_short_ratio=log2((tus$B6_counts+1)/(tus$CAST_counts+1))
tus$log2_long_short_ratio[tus$AT_long=="CAST"] = log2((tus$CAST_counts+1)/(tus$B6_counts+1))[tus$AT_long=="CAST"] 

load(paste(organ, "_gencode.vM25.annotation.gene-Readcounts.RData", sep=""))
View(merge.counts)
dim(merge.counts)
tus_gene = merge(tus , merge.counts, by = "geneID")
dim(tus_gene)

plot(tus_gene$log2_B6_CAST_ratio, log2((tus_gene$B6.proseq+1)/(tus_gene$CAST.proseq+1)),
     xlab="log2((B6+1)/(CAST+1)) Proseq in Tunit gene body",
     ylab="log2((B6+1)/(CAST+1)) Proseq in gene annotation",
     pch=19, col = rgb(0, 0 , 0, alpha=0.3), main=organ
     )
cor.test(tus_gene$log2_B6_CAST_ratio, log2((tus_gene$B6.proseq+1)/(tus_gene$CAST.proseq+1)),)

boxplot(log2((tus_gene$B6.rna+1)/(tus_gene$CAST.rna+1)) ~ tus_gene$AT_long,
        ylim=c(-1.2,1.2),
        outline=F,las=1,
        main=organ)

pdf(paste(organ,"_AT_GeneAnnotation_RNA-seq-exon-read_ratio.pdf", sep=""), width=7, height = 7)
#par(mfrow=c(1,3))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
# gene annotation, doesn't remove the AT window.

boxplot(log2((tus_gene$B6.exon+1)/(tus_gene$CAST.exon+1)) ~ tus_gene$AT_long,
        ylim=c(-1.2,1.2),
        outline=F,las=1,
        main=organ)

abline(h=0, col="red")
dev.off()

sum(tus_gene$AT_long=="B6")
sum(tus_gene$AT_long=="CAST")
sum(tus_gene$group=="NoAT")
temp=log2((tus_gene$B6.exon+1)/(tus_gene$CAST.exon+1))
wilcox.test(tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="B6"], 
            tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="CAST"])   
wilcox.test(tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="B6"], 
            tus$log2_B6_CAST_ratio[tus$group=="NoAT"])  
wilcox.test(tus$log2_B6_CAST_ratio[tus$group=="NoAT"], 
            tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="CAST"])




if(0){
library("vioplot")
vioplot(log2((tus$B6_counts+1)/(tus$CAST_counts+1)) ~ tus$AT_long,
        main=organ)
abline(h=0, col="red")
boxplot(log2((tus$B6_counts+1)/(tus$CAST_counts+1)) ~ tus$AT_long,
        main=organ)
abline(h=0, col="red")

vioplot(tus$log2_long_short_ratio ~ tus$group,
        main=organ)
abline(h=0, col="red")
boxplot(tus$log2_long_short_ratio ~ tus$group,
        main=organ)
abline(h=0, col="red")

cat("NoAT, N=" , sum(tus$group=="NoAT"))
cat("WithAT, N=" , sum(tus$group=="AT"))
}
# withAT
# Tunit need to remove reads in AT window
tus_withAT=read.table(file = paste(organ,"_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed",sep=""), header=FALSE)
colnames(tus_withAT)[4]="intersect"
colnames(tus_withAT)[16]="tunit"
colnames(tus_withAT)[19]="AT_long"
colnames(tus_withAT)[14] = "TXSTART"
colnames(tus_withAT)[15] = "TXEND"
colnames(tus_withAT)[18] = "TXSTRAND"
tus_withAT <- tus_withAT[(tus_withAT$TXEND-tus_withAT$TXSTART)>500,]
dim(tus_withAT)
# remove dupplicated tus, keep one of the intersect
tus_withAT = tus_withAT[!duplicated(tus_withAT$tunit),]
dim(tus_withAT)
b6.bw=paste(organ, "_map2ref_1bpbed_map5_B6", sep="")
cast.bw=paste(organ, "_map2ref_1bpbed_map5_CAST", sep="")
filenames <- c(b6.bw, cast.bw)

bodies <- tus_withAT[,13:18] # tunit
colnames(bodies)[6]="TXSTRAND"
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+250
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-250
## Gets counts
counts <- NULL
for(f in filenames) {
    counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE, path=file_dir))
}


tus_withAT$tunit_B6_counts = counts[,1]
tus_withAT$tunit_CAST_counts = counts[,2]

bodies <- tus_withAT[,1:6]  #intersect region, AT window
## Gets counts
counts <- NULL
for(f in filenames) {
    counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE, path=file_dir))
}

tus_withAT$it_B6_counts = counts[,1]
tus_withAT$it_CAST_counts = counts[,2]

tus_withAT$B6_counts = tus_withAT$tunit_B6_counts - tus_withAT$it_B6_counts
tus_withAT$CAST_counts =tus_withAT$tunit_CAST_counts - tus_withAT$it_CAST_counts
tus_withAT$log2_B6_CAST_ratio=log2((tus_withAT$B6_counts+1)/(tus_withAT$CAST_counts+1))
tus_withAT$log2_long_short_ratio=log2((tus_withAT$B6_counts+1)/(tus_withAT$CAST_counts+1))
tus_withAT$log2_long_short_ratio[tus_withAT$AT_long=="CAST"] = log2((tus_withAT$CAST_counts+1)/(tus_withAT$B6_counts+1))[tus_withAT$AT_long=="CAST"] 

if(0){
hist(tus_withAT$log2_long_short_ratio, 
     breaks = seq(-5.0025,5,0.05), 
     main=organ)
abline(v=0, col="red")

vioplot(tus_withAT$log2_long_short_ratio, tus$log2_long_short_ratio[tus$group=="AT"],tus$log2_long_short_ratio[tus$group=="NoAT"],
        names=c("exclude AT", "include AT", "tunit without AT"),
        ylab="log2(long+1/short+1)",
        main=organ)
abline(h=0, col="red")
}
pdf(paste(organ,"_AT_Proseq_read_ratio.pdf", sep=""), width=7, height = 7)
#par(mfrow=c(1,3))
par(mar=c(6.1, 7.1, 2.1, 2.1)) #d l u r 5.1, 4.1, 4.1, 2.1
par(mgp=c(3,1,0))
par(cex.lab=2.2, cex.axis=2.2)
boxplot(tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="B6"], 
        tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="CAST"],
        tus$log2_B6_CAST_ratio[tus$group=="NoAT"],
        names=c("B6", "CAST", "No AT"),
        xlab="AT long allele",
        ylab="log2(B6+1/CAST+1)",
        ylim=c(-1.2,1.2),
        outline=F,las=1,
        main=organ)
abline(h=0, col="red") # remove the outliner
dev.off()

sum(tus_withAT$AT_long=="B6")
sum(tus_withAT$AT_long=="CAST")
sum(tus$group=="NoAT")
wilcox.test(tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="B6"], 
            tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="CAST"])   
wilcox.test(tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="B6"], 
            tus$log2_B6_CAST_ratio[tus$group=="NoAT"])  
wilcox.test(tus$log2_B6_CAST_ratio[tus$group=="NoAT"], 
            tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="CAST"])



plot(ecdf(tus$log2_B6_CAST_ratio[tus$group=="NoAT"]),
     #pch=10, 
     col="dark blue",
     #xlim=c(0,2),
     #ylim=c(0.6,1),
     xlab="log2(B6+1/CAST+1)",
     ylab="CDF",
     las=1,cex=1,
     main=organ
)
plot(ecdf(tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="B6"]),
     #pch=10, 
     col="red",
     #xlim=c(0,2),
     #ylim=c(0.6,1),
     xlab="log2(B6+1/CAST+1)",
     ylab="CDF",
     las=1,cex=1,
     main=organ
)
lines(ecdf(tus_withAT$log2_B6_CAST_ratio), cex=1, col="dark orange")
lines(ecdf(tus$log2_B6_CAST_ratio[tus$group=="NoAT"]), cex=1, col="darkblue")


lines(ecdf(tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="B6"]), cex=1, col="red")

lines(ecdf(tus_withAT$log2_B6_CAST_ratio[tus_withAT$AT_long=="CAST"]), cex=1, col="dark green")

legend("right", 
       legend = c("B6","CAST" ),
       title = ,
       #pch=c(15,15),
       #lty=c(0,0),
       #bty="n",
       lwd=1.5, 
       #density=c(10000,25),
       #angle=c(180,45),
       #angle=45,
       #fill=c("blue","dark organe","dark green")
       col=c("red","dark green"),
       pch=c(19,19),
       bty = "n"
)




# relationship between pro-seq and rna-seq 
# from F1_RNAseq_getCounts.R
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/RNA-seq")
load("BN_gencode.vM25.annotation.gene-Readcounts.RData")

counts = merge.counts
counts = counts[counts$All.proseq > 10,]
counts = counts[counts$B6.proseq > 10,]
counts = counts[counts$CAST.proseq > 10,]
counts = counts[counts$All.rna > 10,]
counts = counts[counts$B6.rna > 10,]
counts = counts[counts$CAST.rna > 10,]

plot(counts$All.proseq, counts$All.rna)
plot(log(counts$All.proseq), log(counts$All.rna))
cor.test(log(counts$All.proseq), log(counts$All.rna))
plot(counts$B6.proseq, counts$B6.rna)
plot(log(counts$B6.proseq), log(counts$B6.rna))


#plot(counts$All.exon, counts$All.rna)
#cor(counts$All.exon, counts$All.rna)
plot(counts$All.proseq, counts$All.rna)
cor.test(counts$All.proseq, counts$All.rna)
plot(log(counts$All.proseq+1), log(counts$All.rna+1))
cor.test(log(counts$All.proseq+1), log(counts$All.rna+1))
plot(counts$All.proseq, counts$All.exon)
cor.test(counts$All.proseq, counts$All.exon)
plot(log(counts$All.proseq+1), log(counts$All.exon+1))
cor.test(log(counts$All.proseq+1), log(counts$All.exon+1))


