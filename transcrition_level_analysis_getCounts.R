##
## getCounts.R - Counts reads in each gene.
require(bigWig)

tus <- read.table("./tunit_protein_coding/BN_all_h5.preds.full_inProtein_coding.bed", header=FALSE)
colnames(tus)[2] = "TXSTART"
colnames(tus)[3] = "TXEND"
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


