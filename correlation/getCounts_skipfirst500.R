#cd /workdir/sc2457/F1_Tissues/DEseq_analysis
#ln -s /workdir/sc2457/mouse_AlleleSpecific/mouse_genome.sanger.ac.uk/REL-1505-SNPs_Indels/PersonalGenome_P.CAST_M.B6_indelsNsnps_CAST.bam/gencode.vM25.annotation.gtf .
#echo -e "#chr\tchrStart\tchrEnd\ttranscript_id\tgene_name\tstrand" > gencode.vM25.annotation.transcript.bed
#cat gencode.vM25.annotation.gtf | awk 'BEGIN{OFS="\t"} $3=="transcript" {split($12,a,"\""); split($16,b,"\""); print $1, $4-1, $5,a[2], b[2],$7}' >> gencode.vM25.annotation.transcript.bed
#cat gencode.vM25.annotation.gtf | awk 'BEGIN{OFS="\t"} $3=="gene" {split($10,a,"\""); split($14,b,"\""); print $1, $4-1, $5,a[2], b[2],$7}' > gencode.vM25.annotation.gene.bed

## getCounts.R - Counts reads in each gene.
require(bigWig)

tus <- read.table("gencode.vM25.annotation.transcript.bed", header=F)
tus <- tus[(tus$V3-tus$V2)>10000,]  #use TRX size >10K


#geneID_name <- read.table("gencode.vM20_geneID_name_pair.txt", header=F)
#colnames(geneID_name)=c("GENEID", "GENENAME")

# add geneID and gene name to the table
#library("sqldf")
#new_tus <- sqldf ("select V1,V2,V3,GENEID,GENENAME,V6 from tus left join geneID_name on tus.V4=geneID_name.GENEID" )
#bodies <- new_tus
bodies=tus
bodies <- bodies[(bodies$V1 != "chrM"), ]

bodies$V2[bodies$V6 == "+"] <- bodies$V2[bodies$V6 == "+"]+500
bodies$V3[bodies$V6 == "-"] <- bodies$V3[bodies$V6 == "-"]-500
bodies <-unique(bodies)


countBigWig <- function(prefix, bed, rpkm=FALSE, path="./") {
 pl <- load.bigWig(paste(path, prefix, "_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_minus.bw", sep=""))

 counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
        if(rpkm==TRUE) {
                counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }

 return(counts)
}

stage     <- c("BN", "GI", "HT", "LV", "SK", "SP", "ST", "KD")
replicate <- c("MB6_A", "MB6_F", "MB6_G", "PB6_B", "PB6_C", "PB6_D", "PB6_E")

filenames=c()
for (s in stage){
	filenames <- c(filenames, paste(s, replicate, sep="_"))
}

if(0){
stage     <- c("HT", "SK", "KD")
replicate <- c("PB6_F5","PB6_F6")

for (s in stage){
	filenames <- c(filenames, paste(s, replicate, sep="_"))
}
}

## Gets counts
counts <- NULL
for(f in filenames) {
	counts <- cbind(counts, countBigWig(f, bodies, rpkm=FALSE))
}
colnames(counts) <- filenames
save.image("data-counts.vM25.transcript_withSingle.RData") #transcript
#save.image("data-counts.vM25.transcript.RData") #transcript
#save.image("data-counts.vM25.gene.RData") #gene
remove(counts)


## Gets RPKMs
rpkm <- NULL
for(f in filenames) {
    rpkm <- cbind(rpkm, countBigWig(f, bodies, rpkm=TRUE))
}
colnames(rpkm) <- filenames

#save.image("data-rpkms.vM25.transcript.RData")
save.image("data-rpkms.vM25.transcript_withSingle.RData")
#save.image("data-rpkms.vM25.gene.RData")
remove(rpkm)


## make cluster
load("data-counts.vM25.RData")
load("data-rpkms.vM25.RData")
min_count = 5
colnames_keep = colnames(counts)[grep("LG", colnames(counts),invert = T)]
colnames_keep = colnames_keep[grep("TH", colnames_keep,invert = T)]
counts <- counts[,colnames_keep ]
indx_counts <- rowSums(counts>=min_count) >= dim(counts)[2]  #every sample has at least min_count reads
indx_trxSize<- (bodies[,3]-bodies[,2])>10000  # to get a robost signal
indx <- indx_counts & indx_trxSize
sub_rpkm <- rpkm[indx,colnames_keep]

write.table(sub_rpkm , file = "rpkm_5reads_trx10K_bodyafter500bp_noLG.txt", quote =FALSE, sep="\t") #full length

## make cluster WITHOUT single base run-on 
load("data-counts.vM25.transcript.RData")
load("data-rpkms.vM25.transcript.RData")
min_count = 5
colnames_keep = colnames(counts)[grep("LG", colnames(counts),invert = T)]
colnames_keep = colnames_keep[grep("TH", colnames_keep,invert = T)]
colnames_keep = colnames_keep[grep("F5", colnames_keep,invert = T)]
colnames_keep = colnames_keep[grep("F6", colnames_keep,invert = T)]

counts <- counts[,colnames_keep ]
indx_counts <- rowSums(counts>=min_count) >= dim(counts)[2]  #every sample has at least min_count reads
indx_trxSize<- (bodies[,3]-bodies[,2])>10000  # to get a robost signal
indx <- indx_counts #& indx_trxSize
sub_rpkm <- rpkm[indx,colnames_keep]

write.table(sub_rpkm , file = "rpkm_5reads_trx10K_bodyafter500bp_noLG_noSingleBase.vM25.transcript.txt", quote =FALSE, sep="\t", col.names = TRUE) 
write.table(sub_rpkm , file = "rpkm_5reads_trx10K_bodyafter500bp_noLG_noSingleBase.vM25.gene.txt", quote =FALSE, sep="\t", col.names = TRUE) 

## make cluster WITH single base run-on 
load("data-counts.vM25.transcript_withSingle.RData")
load("data-rpkms.vM25.transcript_withSingle.RData")
min_count = 5
colnames_keep = colnames(counts)[grep("LG", colnames(counts),invert = T)]
colnames_keep = colnames_keep[grep("TH", colnames_keep,invert = T)]

counts <- counts[,colnames_keep ]
indx_counts <- rowSums(counts>=min_count) >= dim(counts)[2]  #every sample has at least min_count reads
indx_trxSize<- (bodies[,3]-bodies[,2])>10000  # to get a robost signal
indx <- indx_counts #& indx_trxSize
sub_rpkm <- rpkm[indx,colnames_keep]

write.table(sub_rpkm , file = "rpkm_5reads_trx10K_bodyafter500bp_withSingleBase.vM25.transcript.txt", quote =FALSE, sep="\t", col.names = TRUE) 


