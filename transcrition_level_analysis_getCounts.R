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
#colnames(tus)[4] = "TXNAME"
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
vioplot(log2((tus$B6_counts+1)/(tus$CAST_counts+1)) ~ tus$V12,
        main=organ)
abline(h=0)
boxplot(log2((tus$B6_counts+1)/(tus$CAST_counts+1)) ~ tus$V12,
        main=organ)
abline(h=0)

cat("g9, N=" , sum(tus$V12=="g9"))
cat("m, N=" , sum(tus$V12=="m"))
cat("s, N=" , sum(tus$V12=="s"))

change_ratio=2
testor <-  matrix(c(sum(abs(tus$change)>change_ratio & tus$V12=="m"), sum(tus$V12=="m"),
                    sum(abs(tus$change)>change_ratio & tus$V12=="g9"),sum(tus$V12=="g9")),
                  nrow = 2,
                  dimnames = list(TSS = c("change>1", "all"),
                                  KS.Test = c("m", "g9"))); testor

testor <-  matrix(c(sum(abs(tus$change)>change_ratio & tus$V12=="s"), sum(tus$V12=="s"),
                    sum(abs(tus$change)>change_ratio & tus$V12=="g9"),sum(tus$V12=="g9")),
                  nrow = 2,
                  dimnames = list(TSS = c("change>1", "all"),
                                  KS.Test = c("s", "g9"))); testor
change_ratio=2
testor <-  matrix(c(sum(abs(tus$change)>change_ratio & (tus$V12=="s"|tus$V12=="m")), sum(tus$V12=="s"|tus$V12=="m"),
                    sum(abs(tus$change)>change_ratio & tus$V12=="g9"),sum(tus$V12=="g9")),
                  nrow = 2,
                  dimnames = list(log2.change = c(paste("change>",change_ratio, sep=""), "all"),
                                  KS.Test = c("s+m", "g9"))); testor

f = fisher.test(testor, alternative = "g" ); f
#p.value= c(p.value, f$p.value)
#odds.ratio = c(odds.ratio, f$estimate)


colnames(counts) <- filenames
counts <- cbind.data.frame(counts, bodies$TXNAME)
plot(1/2*(log2(counts$BN_map2ref_1bpbed_map5_B6)+log2(counts$BN_map2ref_1bpbed_map5_CAST)),
     log2(counts$BN_map2ref_1bpbed_map5_B6+1/counts$BN_map2ref_1bpbed_map5_CAST+1))

View(counts)
dim(bodies)



