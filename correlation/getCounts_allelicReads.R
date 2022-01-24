##
## getCounts.R - Counts reads in each gene.
# cd /workdir/sc2457/F1_Tissues/3rd_batch/map2ref_1bpbed_map5
require(bigWig)

tus <- read.table("gencode.vM25.annotation.gene.bed", header=F)
tus <- read.table("gencode.vM25.annotation.gene_inStrainDomain.bed", header=F)
tus <- read.table("gencode.vM25.annotation.gene_inImprintingDomain.bed", header=F)
tus <- tus[(tus$V3-tus$V2)>10000,]  #500 for gene body ##10K to get a robust signal
bodies <- tus
bodies <- bodies[(bodies$V1 != "chrM"), ]
bodies <- bodies[(bodies$V1 != "chrX"), ]
bodies <- bodies[(bodies$V1 != "chrY"), ]



countBigWig_mat <- function(prefix, bed, rpkm=FALSE, path="./") {
 pl <- load.bigWig(paste(path, prefix, "_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_R1.mat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted_minus.bw", sep=""))

 counts <- bed6.region.bpQuery.bigWig(pl, mn, bed, abs.value = TRUE)
        if(rpkm==TRUE) {
                counts <- counts * (1000/(bed[,3]-bed[,2])) * (1e6/(abs(pl$mean)*pl$basesCovered+abs(mn$mean)*mn$basesCovered))
        }

 return(counts)
}

countBigWig_pat <- function(prefix, bed, rpkm=FALSE, path="./") {
 pl <- load.bigWig(paste(path, prefix, "_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted_plus.bw", sep=""))
 mn <- load.bigWig(paste(path, prefix, "_R1.pat.bowtie.gz_AMBremoved_sorted_specific.map2ref.map5.1bp.sorted_minus.bw", sep=""))

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


## Gets counts
mat_counts <- NULL
pat_counts <- NULL
for(f in filenames) {
	mat_counts <- cbind(mat_counts, countBigWig_mat(f, bodies, rpkm=FALSE))
    pat_counts <- cbind(pat_counts, countBigWig_pat(f, bodies, rpkm=FALSE))
}
head(mat_counts)
colnames(mat_counts) <- filenames
colnames(pat_counts) <- filenames
rownames(mat_counts) <- bodies$V4
rownames(pat_counts) <- bodies$V4

## Create DESeq2 object.
library("DESeq2")
organ <- c(rep("BN", 7), rep("GI", 7), rep("HT", 7), rep("LV", 7), 
           rep("SK", 7), rep("SP", 7), rep("ST", 7), rep("KD", 7))
condition <- c(rep("MB6", 3), rep("PB6", 4))
replicate <- factor(rep(c(1:3, 1:4), 8))
Design <- data.frame(colnames(mat_counts), organ, condition, replicate)
dds <- DESeqDataSetFromMatrix(countData= mat_counts, colData= Design, design= ~organ+condition+replicate)
vsd <- vst(dds, blind=FALSE)
# for imprinted due to small n
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
#
mat_counts_vsd = assay(vsd)

dds_p <- DESeqDataSetFromMatrix(countData= pat_counts, colData= Design, design= ~organ+condition+replicate)
vsd_p <- vst(dds_p, blind=FALSE)
# for imprinted due to small n
vsd_p <- varianceStabilizingTransformation(dds_p, blind=FALSE)
#
pat_counts_vsd = assay(vsd_p)

save.image("data-counts.RData")
save.image("data-counts-strain.RData")
#save.image("data-counts-imprinting.RData")
#remove(counts)

## make cluster
setwd("/Users/sc2457/Box Sync/Danko_lab_work/F1_8Tissues/correlation/allelicReads/")
load("data-counts.RData")
load("data-counts-strain.RData")
load("data-counts-imprinting.RData")
#View(mat_counts)
#View(mat_counts_vsd)
#View(pat_counts)
#counts = mat_counts + pat_counts
View(counts)
if(0){
min_count = 1
counts = mat_counts + pat_counts
indx_counts <- rowSums(counts>=min_count) >= dim(counts)[2]  #every sample has at least min_count reads
indx_trxSize<- (bodies[,3]-bodies[,2])>10000  # to get a robost signal
indx <- indx_counts & indx_trxSize
mat_counts = mat_counts[indx,]
pat_counts = pat_counts[indx,]
}

p_value = pat_counts

for (i in (1:dim(mat_counts)[1])){
    for (j in (1:dim(mat_counts)[2])){
        x=mat_counts[i,j]
        n=mat_counts[i,j] + pat_counts[i,j]
        if (n>0){b=binom.test(x, n, p = 0.5,
                              alternative = c("two.sided"),
                              conf.level = 0.95)
        p_value[i,j]=b$p.value} else {
            p_value[i,j]=1
        }
        
    }
}
#View(p_value)
dim(p_value)
fdr <- p_value

for (j in (1:dim(mat_counts)[2])){
    fdr[,j] = p.adjust(p_value[,j], method = "fdr", n = length(p_value[,j]))
}
#View(fdr)

fdr_seven <- fdr
for (i in (1:dim(fdr)[1])){
    for (j in seq(1,dim(fdr)[2],7)){
        fdr_seven[i,j:(j+6)] = max(fdr[i,j:(j+6)])
    }
}
#View(fdr_seven)
if(0){
    library("metap")
    pvalue_fs <- p_value
    for (i in (1:dim(p_value)[1])){
        for (j in seq(1,dim(p_value)[2],7)){
            pvalue_fs[i,j:(j+6)] = sumlog(p_value[i,j:(j+6)])$p
        }
    }
    View(pvalue_fs)
}

# only keep genes with allelic biased in all replicates of a organ
f = apply(fdr_seven, 1, min)
f = apply(fdr, 1, min)
sum(f<0.1)
#mat_counts = mat_counts[which(f<0.1),]
#pat_counts = pat_counts[which(f<0.1),]

new_mat_counts_vsd = mat_counts_vsd[which(f<0.1),]
new_pat_counts_vsd = pat_counts_vsd[which(f<0.1),]
allelic_diff = new_mat_counts_vsd - new_pat_counts_vsd

new_mat_counts = mat_counts[which(f<0.1),]
new_pat_counts = pat_counts[which(f<0.1),]
allelic_diff = log2((new_mat_counts+1)/(new_pat_counts+1))


write.table(allelic_diff , file = "allelic_diff_trx10K_vsd_fdr0.1.txt", quote =FALSE, sep="\t") #full length
write.table(allelic_diff , file = "allelic_diff_trx10K_vsd_fdrseven0.1.txt", quote =FALSE, sep="\t") #full length

write.table(allelic_diff , file = "allelic_diff_trx10K_strain_vsd_fdr0.1.txt", quote =FALSE, sep="\t") #full length
write.table(allelic_diff , file = "allelic_diff_trx10K_strain_vsd_fdrseven0.1.txt", quote =FALSE, sep="\t") #full length

write.table(allelic_diff , file = "allelic_diff_trx10K_imprinting_vsd_fdr0.1.txt", quote =FALSE, sep="\t") #full length
write.table(allelic_diff , file = "allelic_diff_trx10K_imprinting_vsd_fdrseven0.1.txt", quote =FALSE, sep="\t") #full length

write.table(allelic_diff , file = "allelic_diff_trx10K_fdrseven0.1.txt", quote =FALSE, sep="\t") #full length
write.table(allelic_diff , file = "allelic_diff_trx10K_fdr0.1.txt", quote =FALSE, sep="\t") #full length

write.table(allelic_diff , file = "allelic_diff_trx10K_strain_fdr0.1.txt", quote =FALSE, sep="\t") #full length
write.table(allelic_diff , file = "allelic_diff_trx10K_strain_fdrseven0.1.txt", quote =FALSE, sep="\t") #full length

write.table(allelic_diff , file = "allelic_diff_trx10K_imprinting_fdr0.1.txt", quote =FALSE, sep="\t") #full length
write.table(allelic_diff , file = "allelic_diff_trx10K_imprinting_fdrseven0.1.txt", quote =FALSE, sep="\t") #full length

#write.table(allelic_diff , file = "allelic_diff_1reads_trx10K.txt", quote =FALSE, sep="\t") #full length
#write.table(allelic_diff , file = "allelic_diff_1reads_trx10K_strain.txt", quote =FALSE, sep="\t") #full length
#write.table(allelic_diff , file = "allelic_diff_1reads_trx10K_imprinting.txt", quote =FALSE, sep="\t") #full length
#write.table(allelic_diff , file = "allelic_diff_trx10K_strain_fdr7.txt", quote =FALSE, sep="\t") #full length
#write.table(allelic_diff , file = "allelic_diff_trx10K_strain_pvalue_fs0.05.txt", quote =FALSE, sep="\t") #full length
#write.table(allelic_diff , file = "allelic_diff_trx10K_strain_vsd_fdr7.txt", quote =FALSE, sep="\t") #full length
#write.table(allelic_diff , file = "allelic_diff_trx10K_vsd_fdr7.txt", quote =FALSE, sep="\t") #full length



