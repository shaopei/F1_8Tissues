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
organ="SP"
tus=read.table("BN_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed") 
tus_noAT=read.table(file = paste(organ,"_all_h5.preds.full_inProtein_coding_withoutATwindow.bed",sep=""), header=FALSE)
tus_noAT = tus_noAT[,1:6] 
tus_noAT$AT_long = "NoAT"
tus_noAT$group = "NoAT"
dim(tus_noAT)

tus_withAT=read.table(file = paste(organ,"_AT_4tunitIntersectNativeHMM_intersectRegion_strain.bed",sep=""), header=FALSE)
tus_withAT = tus_withAT[,13:19]  #tunits
tus_withAT$group = "AT"
colnames(tus_withAT)[1:8] = colnames(tus_noAT)[1:8]
tus=rbind(tus_noAT, tus_withAT)
dim(tus)
#View(tus)
tus=unique(tus)

colnames(tus)[2] = "TXSTART"
colnames(tus)[3] = "TXEND"
colnames(tus)[6] = "TXSTRAND"
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

# withAT
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

bodies <- tus_withAT[,13:18]
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

hist(tus_withAT$log2_long_short_ratio, 
     breaks = seq(-5.0025,5,0.05), 
     main=organ)
abline(v=0, col="red")

vioplot(tus_withAT$log2_long_short_ratio, tus$log2_long_short_ratio[tus$group=="AT"],tus$log2_long_short_ratio[tus$group=="NoAT"],
        names=c("exclude AT", "include AT", "tunit without AT"),
        ylab="log2(long+1/short+1)",
        main=organ)
abline(h=0, col="red")
