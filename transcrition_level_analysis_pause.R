file_dir="~/Box Sync/KD_IGV/2020July/"
setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/transcription_level_analysis/pause/")

organ="SK"
tus=read.table(file = paste(organ,"_dREG_5mat5pat_uniq_AllReads_maxPauseinProteinCodingTunit_mat_patSeq.bed",sep=""), header=FALSE)

colnames(tus)[13:14] = c("B6_allele", "CAST_allele")
# keep tus with C in one allele
dim(tus)
tus=tus[(tus$B6_allele == "C" | tus$CAST_allele == "C"),]
dim(tus)
# two group
# One C, one nonC
tus$group = "2C"
tus$nonC=tus$CAST_allele
tus$nonC[tus$CAST_allele=="C"] = tus$B6_allele[tus$CAST_allele=="C"]
tus$group[(tus$nonC!= "C")] = "1C"
dim(tus)
# per 1C/2C group, tunit can only be count once
tus=tus[!duplicated(tus[,c(7,8,9,12,15)]),] 
dim(tus)
sum(tus$group == "1C") 
sum(tus$group == "2C")

#View(tus)

colnames(tus)[8] = "TXSTART"
colnames(tus)[9] = "TXEND"
colnames(tus)[12] = "TXSTRAND"
tus <- tus[(tus$TXEND-tus$TXSTART)>=500,]
dim(tus)
sum(tus$group == "1C") 
sum(tus$group == "2C")
# body exclude the first h bp of transcript
h=300
bodies <- tus[,7:12]
bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+h
bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-h

countBigWig <- function(prefix, bed, rpkm=FALSE, path="./") {
  pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
  mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))
  
  counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
  if(rpkm==TRUE) {
    counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
  }
  
  return(counts)
}

b6.bw=paste(organ, "_PB6_F5N6_map2ref_1bpbed_map3_B6", sep="")
cast.bw=paste(organ, "_PB6_F5N6_map2ref_1bpbed_map3_CAST", sep="")
filenames <- c(b6.bw, cast.bw)
## Gets counts
counts <- NULL
for(f in filenames) {
  counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE, path=file_dir))
}

tus$B6_counts = counts[,1]
tus$CAST_counts = counts[,2]
tus$log2_C_NonC_ratio = log2((tus$B6_counts+1)/(tus$CAST_counts+1))
tus$log2_C_NonC_ratio[tus$CAST_allele == "C"] = log2((tus$CAST_counts+1)/(tus$B6_counts+1))[tus$CAST_allele == "C"]

library("vioplot")
vioplot(tus$log2_C_NonC_ratio ~ tus$group,
        main=paste(organ, " Tunit Proseq", sep=""))
abline(h=0, col="red")
boxplot(tus$log2_C_NonC_ratio ~ tus$group,
        #ylim=c(-1,1),
        #xlab="AT long allele",
        ylab="log2(C+1/NonC+1)",
        outline=F,las=1,
        main=paste(organ, " Tunit Proseq", sep=""))
abline(h=0, col="red")
dim(tus)
sum(tus$group == "1C") 
sum(tus$group == "2C")

wilcox.test(tus$log2_C_NonC_ratio[tus$group=="1C"],
            tus$log2_C_NonC_ratio[tus$group=="2C"])
