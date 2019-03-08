##
## getCounts.R - Counts reads in each gene.
require(bigWig)

tus <- read.table("gencode.vM20.annotation_transcript.bed", header=F)
tus <- tus[(tus$V3-tus$V2)>500,]


geneID_name <- read.table("gencode.vM20_geneID_name_pair.txt", header=F)
colnames(geneID_name)=c("GENEID", "GENENAME")

library("sqldf")
new_tus <- sqldf ("select V1,V2,V3,GENEID,GENENAME,V6 from tus left join geneID_name on tus.V4=geneID_name.GENEID" )
bodies <- new_tus
bodies <- bodies[(bodies$V1 != "chrM"), ]

countBigWig <- function(prefix, bed, rpkm=FALSE, path="./") {
 pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))

 counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
        if(rpkm==TRUE) {
                counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }

 return(counts)
}

stage     <- c("BN", "GI", "HT", "LG", "LV", "SK", "SP", "ST")
replicate <- c("MB6_A", "MB6_F", "MB6_G", "PB6_B", "PB6_C", "PB6_D", "PB6_E")

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
save.image("data-counts.RData")
remove(counts)


## Gets RPKMs
rpkm <- NULL
for(f in filenames) {
    rpkm <- cbind(rpkm, countBigWig(f, bodies, rpkm=TRUE))
}
colnames(rpkm) <- filenames

save.image("data-rpkms.RData")
remove(rpkm)



