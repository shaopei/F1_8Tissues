##
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

filenames <- c("BN_map2ref_1bpbed_map5", "BN_map2ref_1bpbed_map5_B6", "BN_map2ref_1bpbed_map5_CAST")

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

exon.filenames <- c("BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.exon.read.count",
	"BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.nsorted.exon.read.count", 
	"BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast.nsorted.exon.read.count")
exon.counts <- NULL
for(f in exon.filenames) {
	r <- read.table(paste("/workdir/sc2457/F1_Tissues/RNA-seq/STAR_BN/",f,sep=""), header=F)
	#merge.counts <- merge(rna.counts, r, by="V1", all=TRUE) 
	exon.counts <- cbind(exon.counts, r[,2])
}
exon.counts <- cbind.data.frame(r$V1, exon.counts)
colnames(exon.counts) <- c("geneID","All.exon", "B6.exon", "CAST.exon")


rna.filenames <- c("BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.gene.read.count",
	"BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.b6.nsorted.gene.read.count", 
	"BN_MB6_BOTH_RNA_mat3waspAligned.sortedByCoord.out.cast.nsorted.gene.read.count")
rna.counts <- NULL
for(f in rna.filenames) {
	r <- read.table(paste("/workdir/sc2457/F1_Tissues/RNA-seq/STAR_BN/",f,sep=""), header=F)
	#merge.counts <- merge(rna.counts, r, by="V1", all=TRUE) 
	rna.counts <- cbind(rna.counts, r[,2])
}
rna.counts <- cbind.data.frame(r$V1, rna.counts)
colnames(rna.counts) <- c("geneID","All.rna", "B6.rna", "CAST.rna")


merge.counts <- merge(counts, rna.counts, by="geneID", all.x=TRUE) 
merge.counts <- merge(merge.counts, exon.counts, by="geneID", all.x=TRUE) 

save.image("BN_gencode.vM25.annotation.gene-Readcounts.RData")

setwd("~/Box Sync/Danko_lab_work/F1_8Tissues/PolyA_Allele-specific/RNA-seq")
load("BN_gencode.vM25.annotation.gene-Readcounts.RData")

counts = merge.counts
counts$sB6 <- counts$B6.exon/counts$B6.proseq
counts$sCAST <- counts$CAST.exon/counts$CAST.proseq

plot(counts$sB6, counts$sCAST,
     xlim=c(0,200), ylim = c(0,200))






library("sqldf")
new_tus <- sqldf ("select TXCHROM,TXSTART,TXEND,GENEID,TXNAME,TXSTRAND,V2,TXTYPE from tus left join geneID_name on tus.GENEID=geneID_name.V1" )
colnames(new_tus)[7]="GENENAME"
tus <- new_tus
bodies <- tus



tus <- read.table("../annotations/tuSelecter/final_tus.txt", header=TRUE)
tus <- tus[(tus$TXEND-tus$TXSTART)>500,]
geneID_name <- read.table("../annotations/gencode.vM20_geneID_name_pair.txt", header=F)
#colnames(geneID_name)=c("GENEID", "GENENAME")

library("sqldf")
new_tus <- sqldf ("select TXCHROM,TXSTART,TXEND,GENEID,TXNAME,TXSTRAND,V2,TXTYPE from tus left join geneID_name on tus.GENEID=geneID_name.V1" )
colnames(new_tus)[7]="GENENAME"
tus <- new_tus
bodies <- tus
# bodies$TXSTART[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+250
# bodies$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-250

# pause <- tus
# pause$TXEND[bodies$TXSTRAND == "+"] <-bodies$TXSTART[bodies$TXSTRAND == "+"]+250
# pause$TXSTART[bodies$TXSTRAND == "-"] <- bodies$TXEND[bodies$TXSTRAND == "-"]-250

# postcps <- tus
# postcps$TXSTART[bodies$TXSTRAND == "+"] <- bodies$TXEND[bodies$TXSTRAND == "+"]
# postcps$TXEND[bodies$TXSTRAND == "+"] <- bodies$TXEND[bodies$TXSTRAND == "+"]+15000
# postcps$TXEND[bodies$TXSTRAND == "-"] <- bodies$TXSTART[bodies$TXSTRAND == "-"]
# postcps$TXSTART[bodies$TXSTRAND == "-"] <- bodies$TXSTART[bodies$TXSTRAND == "-"]-15000
# postcps$TXSTART[postcps$TXSTART < 0] <- 0

countBigWig <- function(prefix, bed, rpkm=FALSE, path="./") {
 pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))

 counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
        if(rpkm==TRUE) {
                counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }

 return(counts)
}

stage     <- c("WT", "MUT")
replicate <- c(1, 2, 3, 4)

filenames <- c(paste("WT_R", replicate, sep=""), paste("MUT_R", replicate, sep=""))
stage     <- c("WT", "MUT", "PHDHET")
filenames <- c(filenames, "PHDHET_R1", "PHDHET_R2")

## Gets counts
counts <- NULL
#pause_counts <- NULL
#postcps_counts <- NULL
for(f in filenames) {
	counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE))
    #pause_counts <- cbind(pause_counts, countBigWig(f, pause, rpkm=FALSE))
	#postcps_counts <- cbind(postcps_counts, countBigWig(f, postcps, rpkm=FALSE)) 
}
colnames(counts) <- filenames
counts_wPHDHET <-counts
#remove(counts)
save.image("data-counts_wPHDHET.RData")
#save.image("data-counts.RData")
#remove(counts_wPHDHET); #remove(pause_counts); remove(postcps_counts)

## Gets RPKMs
rpkm <- NULL
#pause_rpkm <- NULL
#postcps_rpkm <- NULL
for(f in filenames) {
    rpkm <- cbind(rpkm, countBigWig(f, bodies, rpkm=TRUE))
	#pause_rpkm <- cbind(pause_rpkm, countBigWig(f, pause, rpkm=TRUE))
	#postcps_rpkm <- cbind(postcps_rpkm, countBigWig(f, postcps, rpkm=TRUE))
}
colnames(rpkm) <- filenames
rpkm_wPHDHET <-rpkm
#remove(rpkm)
#save.image("data-rpkms.RData")
save.image("data-rpkms_wPHDHET.RData")
#remove(rpkm_wPHDHET)


